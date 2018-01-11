/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICE6D_MB_MODES_TPP_
#define SRC_SERVICE_GPUENERGYSERVICE6D_MB_MODES_TPP_

#ifdef CUDA

#include <nvToolsExt.h>
#include "GPUEnergyService6D_MB_Modes.h"

#include <cassert>
#include "WorkerBuffer.h"
#include "DataManager.h"
#include "Allocator.h"
#include "RingArray.h"
#include "DataItem.h"
#include "DeviceItem.h"
#include "WorkItem.h"
#include "DeviceProtein.h"
#include "DeviceGridUnion.h"
#include "DeviceParamTable.h"
#include "SimParam.h"

#include "transform_MB_modes.h"
#include "interpolation.h"
#include "interpolation_MB_modes.h"
#include "neighborlist.h"
#include "reduction_MB_modes.h"
#include <algorithm>
#include "macros.h"

// ToDo: remove
#include "Types_6D_Config.h"


namespace as {

template<typename REAL>
GPUEnergyService6D_MB_Modes<REAL>::GPUEnergyService6D_MB_Modes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds) :
	GPUEnergyService<Types_6D_MB_Modes<REAL>>(dataMng), _workerId(0), _deviceIds(deviceIds)
{}

template<typename REAL>
struct GPUEnergyService6D_MB_Modes<REAL>::StageResource {
private:
	using workItem_t = typename GPUEnergyService6D_MB_Modes<REAL>::workItem_t;
public:
	d_GridUnion<REAL> gridRec;
	d_GridUnion<REAL> gridLig[LIGANDS_MAX_NUMBER];
	d_Protein<REAL>* rec;
	d_Protein<REAL>* lig[LIGANDS_MAX_NUMBER];
	d_ParamTable<REAL>* table;
	SimParam<REAL>* simParam;
	workItem_t* item;
};

template<typename REAL>
auto GPUEnergyService6D_MB_Modes<REAL>::createStageResource(workItem_t* item, unsigned const& deviceId) -> StageResource {
	/* item pointers */
//			const auto dofs = item->inputBuffer();
	const auto common = item->common();
//			auto results = item->resultBuffer();

	/* get DataItem pointers */
	StageResource stageResource;

	auto gridRec = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridIdRec, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
	assert(gridRec != nullptr);
	stageResource.gridRec 	= gridRec->getDesc();

	for(int i = 0; i < Common_MB_Modes::numLigands; i++){
		auto gridLig = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridIdLig[i], deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(gridLig != nullptr);

		auto lig = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->ligId[i], deviceId)).get();
		assert(lig != nullptr);

		stageResource.gridLig[i] 	= gridLig->getDesc();
		stageResource.lig[i] 		= &lig->desc;
	}
	auto rec = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->recId, deviceId)).get();
	assert(rec != nullptr);
	stageResource.rec 		= &rec->desc;

	auto table = std::dynamic_pointer_cast<DeviceParamTable<REAL>>(this->_dataMng->get(common->tableId, deviceId)).get();
	assert(table != nullptr);
	stageResource.table 	= &table->desc;

	auto simParam = std::dynamic_pointer_cast<SimParam<REAL>>(this->_dataMng->get(common->paramsId)).get();
	assert(simParam != nullptr);
	stageResource.simParam = simParam;

	stageResource.item 	= item;

	return stageResource;
}

template<typename REAL>
class GPUEnergyService6D_MB_Modes<REAL>::Private {
	using dof_t = typename GPUEnergyService6D_MB_Modes<REAL>::input_t;
	using workItem_t = typename GPUEnergyService6D_MB_Modes<REAL>::workItem_t;

public:

	Private() : stagesMngt(numStages), pipeIdx{0,1}, pipeIdxDof{0,1,2}, numItemsInPipe(0) {
		for (unsigned i = 0; i<numStages; ++i) {
			predicates[0][i] = false;
			predicates[1][i] = false;
		}

		for (int i = 0; i<4; ++i) {
			CUDA_CHECK(cudaStreamCreate(&streams[i]));
		}

		for (int i = 0; i<7; ++i) {
			CUDA_CHECK(cudaEventCreate(&events[i], cudaEventDisableTiming));
		}
	}
	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBufferLig( size_t const& numDOFs, size_t const& numAtoms, size_t const& numLigands, size_t const& ligIdx) {
		const size_t atomBufferSize = numDOFs * numAtoms;
		d_defoLig[ligIdx] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSize));
		for (int i = 0; i < 2; ++i) {
			d_potLig[ligIdx][i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
		}
		d_trafoLig[ligIdx] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3, atomBufferSize, atomBufferSize * numLigands));
	}

	void allocateBufferRec( size_t const& numDOFs, size_t const& numAtoms, size_t const& numLigands) {
			const size_t atomBufferSize = numDOFs*numAtoms;
			for (int i = 0; i < 2; ++i) {
				d_potRec[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
			}
			d_trafoRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3, atomBufferSize, atomBufferSize * numLigands));
		d_defoRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSize));
	}

	void allocateBufferDOF( size_t const& numDOFs, size_t const& dofSize, size_t const& numLigands) {
		for (int i = 0; i < 3; ++i) {
			d_dof[i]    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
		}
		for (int i = 0; i < 2; ++i) {
			d_res[i]    = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(numLigands,dofSize * numDOFs));
			h_res[i]    = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(numLigands,dofSize * numDOFs));
		}
	}


	size_t bufferSizeLig(int i) const {
		return d_defoLig[i].bufferSize();
	}

	size_t bufferSizeRec() const {
		return d_defoRec.bufferSize();
	}

	size_t bufferSizeDOF() const {
		return std::min(d_dof[0].bufferSize(),d_dof[1].bufferSize());
	}

	void addItemAndLigandSize(StageResource const& resc) {
		++numItemsInPipe;
		predicates[pipeIdx[0]][0] = true;
		stagesMngt.push(resc);
	}

	void resizeLigBuffersIfRequired(size_t const& numDOFs, size_t const& numAtomsLig,size_t const& numLigands, size_t const& ligIdx) {
		if (numDOFs * numAtomsLig > bufferSizeLig(ligIdx)) {
			allocateBufferLig(numDOFs, numAtomsLig, numLigands,ligIdx);
		}
	}
	void resizeRecBuffersIfRequired(size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numLigands) {
		if (numDOFs * numAtomsRec > bufferSizeRec()) {
			allocateBufferRec(numDOFs, numAtomsRec, numLigands);
		}
	}
	void resizeDOFBuffersIfRequired(size_t const& numDOFs, size_t const& dofSize, size_t const& numLigands) {
		if (numDOFs > bufferSizeDOF()) {
			allocateBufferDOF(numDOFs, dofSize, numLigands);
		}
	}

	void swapBuffers() {
		std::swap(pipeIdx[0], pipeIdx[1]);
		int 	  tmp = pipeIdxDof[2];
		pipeIdxDof[2] = pipeIdxDof[1];
		pipeIdxDof[1] = pipeIdxDof[0];
		pipeIdxDof[0] = tmp;
	}

	bool pipelineEmpty() const {
		return numItemsInPipe == 0;
	}

	void signalItemPassedLastStage() {
		--numItemsInPipe;
	}

	void resetPrediacatesForIteration() {
		for (unsigned i = 0; i < numStages; ++i) {
			predicates[pipeIdx[0]][i] = false;
		}
	}

	void S0_copyH2D() {
		/* check if new item enters the pipeline */
		if (predicates[pipeIdx[0]][0])
		{
			constexpr unsigned stageId = 0;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

//			//DEBUG
//			for (size_t i = 0; i < it->size(); ++i) {
//				std::cout << it->inputBuffer()[i] << std::endl;
//			}

			cudaVerify(cudaStreamWaitEvent(streams[0], events[2], 0));
			cudaVerify(cudaMemcpyAsync(d_dof[pipeIdxDof[0]].get(0), it->inputBuffer(),
						it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, streams[0]));
			cudaVerify(cudaEventRecord(events[0], streams[0]));

		}
	}

	void configureDevice() {

		if (blockSizeReduce == 0) {
			int id;
			cudaVerify(cudaGetDevice(&id));
			cudaDeviceProp deviceProp;
			cudaVerify(cudaGetDeviceProperties(&deviceProp, id));
			size_t sharedMem = deviceProp.sharedMemPerBlock;
			size_t pow2 = 16;
			while (pow2* (Common_Modes::numModesLig + Common_Modes::numModesRec + 13)*sizeof(REAL) < sharedMem) {
				pow2 *= 2;
			}

			blockSizeReduce = pow2 / 2;
		}

	}

	void S1_transform_potForce() {
		/* check if stage 0 was executed in last iteration */
		if (predicates[pipeIdx[1]][0])
		{

			constexpr unsigned stageId = 1;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;


			const unsigned numElRec = it->size() * stageResc.rec->numAtoms;
			assert(numElRec <= bufferSizeRec());
			size_t gridSizeRec = ( numElRec + BLSZ_TRAFO - 1) / BLSZ_TRAFO;

			/* Device: Wait for completion of copyH2D of DOFs to complete */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[0], 0));

			d_DOFPos_receptor(
				BLSZ_TRAFO,
				gridSizeRec,
				streams[2],
				Common_MB_Modes::numLigands,
				stageResc.rec->xPos,
				stageResc.rec->yPos,
				stageResc.rec->zPos,
				stageResc.rec->xModes,
				stageResc.rec->yModes,
				stageResc.rec->zModes,
				d_dof[pipeIdxDof[1]].get(0),
				stageResc.rec->numAtoms,
				stageResc.rec->numModes,
				it->size(),
				d_defoRec.getX(),
				d_defoRec.getY(),
				d_defoRec.getZ(),
				d_trafoRec.getX(),
				d_trafoRec.getY(),
				d_trafoRec.getZ()
				);

			for(unsigned int lig = 0; lig < Common_MB_Modes::numLigands; lig++){

				const unsigned numElLig = it->size() * stageResc.lig[lig]->numAtoms;
				assert(numElLig <= bufferSizeLig(lig));

				size_t gridSizeLig = ( numElLig + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
				/* Perform cuda kernel calls */

				d_DOFPos(
					BLSZ_TRAFO,
					gridSizeLig,
					streams[2],
					Common_MB_Modes::numLigands,
					lig,
					stageResc.lig[lig]->xPos,
					stageResc.lig[lig]->yPos,
					stageResc.lig[lig]->zPos,
					stageResc.lig[lig]->xModes,
					stageResc.lig[lig]->yModes,
					stageResc.lig[lig]->zModes,
					d_dof[pipeIdxDof[1]].get(0),
					stageResc.lig[lig]->numAtoms,
					stageResc.lig[lig]->numModes,
					it->size(),
					d_defoLig[lig].getX(),
					d_defoLig[lig].getY(),
					d_defoLig[lig].getZ(),
					d_trafoLig[lig].getX(),
					d_trafoLig[lig].getY(),
					d_trafoLig[lig].getZ()
					);

			}
				// DEBUG
	//			cudaDeviceSynchronize();
	//			size_t bufferSize = d_trafoLig.bufferSize();
	//			WorkerBuffer<REAL> h_trafoLig(3,bufferSize);
	//			size_t cpySize = h_trafoLig.bufferSize()*sizeof(REAL);
	//
	//			std::cout << "bufferSize: " << bufferSize << " cpySize: " << cpySize << std::endl;
	//			cudaMemcpy(h_trafoLig.getX(),d_trafoLig.getX(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_trafoLig.getY(),d_trafoLig.getY(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_trafoLig.getZ(),d_trafoLig.getZ(), cpySize, cudaMemcpyDeviceToHost);
	//			for(size_t i = 0; i < bufferSize; ++i) {
	//				std::cout << h_trafoLig.getX()[i] << " " << h_trafoLig.getY()[i] << " " << h_trafoLig.getZ()[i] << std::endl;
	//			}
	//			exit(EXIT_SUCCESS);


				//ToDo how to deal with this in mb calculation:
				/* Device: Signal event when transformation has completed */
				cudaVerify(cudaEventRecord(events[2], streams[2]));
				/* Device: Wait for completion of reduction of the previous round */
				cudaVerify(cudaStreamWaitEvent(streams[2], events[5+pipeIdx[1]], 0));



			for(unsigned int lig = 0; lig < Common_MB_Modes::numLigands; lig++){
				const unsigned numElLig = it->size()*stageResc.lig[lig]->numAtoms;
				assert(numElLig <= bufferSizeLig(lig));

				/* Perform cuda kernel calls */
				size_t gridSizeLig = ( numElLig + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
				/* Perform cuda kernel calls */


				d_potForce (
					BLSZ_INTRPL,
					gridSizeLig,
					streams[2],
					stageResc.gridRec.inner,
					stageResc.gridRec.outer,
					*stageResc.lig[lig],
					it->size(),
					d_trafoLig[lig].getX(lig),
					d_trafoLig[lig].getY(lig),
					d_trafoLig[lig].getZ(lig),
					d_potLig[lig][pipeIdx[1]].getX(),
					d_potLig[lig][pipeIdx[1]].getY(),
					d_potLig[lig][pipeIdx[1]].getZ(),
					d_potLig[lig][pipeIdx[1]].getW()); // OK

				d_potForce (
					BLSZ_INTRPL,
					gridSizeRec,
					streams[2],
					d_dof[pipeIdxDof[1]].get(0),
					lig,
					stageResc.gridLig[lig].inner,
					stageResc.gridLig[lig].outer,
					*stageResc.rec,
					it->size(),
					d_trafoRec.getX(lig),
					d_trafoRec.getY(lig),
					d_trafoRec.getZ(lig),
					d_potRec[pipeIdx[1]].getX(),
					d_potRec[pipeIdx[1]].getY(),
					d_potRec[pipeIdx[1]].getZ(),
					d_potRec[pipeIdx[1]].getW()); // OK



				for(unsigned int liglig = 0; liglig < Common_MB_Modes::numLigands; liglig++){
					if(liglig != lig){
						d_potForce (
							BLSZ_INTRPL,
							gridSizeLig,
							streams[2],
							d_dof[pipeIdxDof[1]].get(0),
							liglig,
							stageResc.gridLig[liglig].inner,
							stageResc.gridLig[liglig].outer,
							*stageResc.lig[lig],
							it->size(),
							d_trafoLig[lig].getX(liglig),
							d_trafoLig[lig].getY(liglig),
							d_trafoLig[lig].getZ(liglig),
							d_potLig[lig][pipeIdx[1]].getX(),
							d_potLig[lig][pipeIdx[1]].getY(),
							d_potLig[lig][pipeIdx[1]].getZ(),
							d_potLig[lig][pipeIdx[1]].getW()); // OK


						//ToDO: check if rotation of forces is right. so far the force is only rotated back by the angle of liglig
						}
				}
				// Debug
	//			cudaDeviceSynchronize();
	////			exit(EXIT_SUCCESS);
	//			WorkerBuffer<REAL> h_potLig(4,stageResc.lig->numAtoms);
	//			size_t cpySize = stageResc.lig->numAtoms*sizeof(REAL);
	//			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[1]].getX(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[1]].getY(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[1]].getW(), cpySize, cudaMemcpyDeviceToHost);
	//			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
	////			for(size_t i = 0; i < 20; ++i) {
	//				std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;
	//			}
	//			exit(EXIT_SUCCESS);


	//			IP.d_NLPotForce<false>(it->devLocGridId(), it->devLocRecId(), it->devLocLigId(),it->size(),
	//					&d_trafoLig, d_potLig[pipeIdx[1]],streams[2]);

				d_NLPotForce(
					BLSZ_INTRPL,
					gridSizeLig,
					streams[2],
					stageResc.gridRec.NL,
					*stageResc.rec,
					*stageResc.lig[lig],
					*stageResc.table,
					*stageResc.simParam,
					it->size(),
					d_trafoLig[lig].getX(lig),
					d_trafoLig[lig].getY(lig),
					d_trafoLig[lig].getZ(lig),
					d_potLig[lig][pipeIdx[1]].getX(),
					d_potLig[lig][pipeIdx[1]].getY(),
					d_potLig[lig][pipeIdx[1]].getZ(),
					d_potLig[lig][pipeIdx[1]].getW()); // OK

				//get nl forces on receptor
				d_NLPotForce(
					BLSZ_INTRPL,
					gridSizeRec,
					streams[2],
					stageResc.gridLig[lig].NL,
					*stageResc.lig[lig],
					*stageResc.rec,
					*stageResc.table,
					*stageResc.simParam,
					it->size(),
					d_trafoRec.getX(lig),
					d_trafoRec.getY(lig),
					d_trafoRec.getZ(lig),
					d_potRec[pipeIdx[1]].getX(),
					d_potRec[pipeIdx[1]].getY(),
					d_potRec[pipeIdx[1]].getZ(),
					d_potRec[pipeIdx[1]].getW()
					); // OK
				for(unsigned int liglig = 0; liglig < Common_MB_Modes::numLigands; liglig++){
					if(liglig != lig){
						d_NLPotForce(
							BLSZ_INTRPL,
							gridSizeLig,
							streams[2],
							stageResc.gridLig[liglig].NL,
							*stageResc.rec,
							*stageResc.lig[lig],
							*stageResc.table,
							*stageResc.simParam,
							it->size(),
							d_trafoLig[lig].getX(liglig),
							d_trafoLig[lig].getY(liglig),
							d_trafoLig[lig].getZ(liglig),
							d_potLig[lig][pipeIdx[1]].getX(),
							d_potLig[lig][pipeIdx[1]].getY(),
							d_potLig[lig][pipeIdx[1]].getZ(),
							d_potLig[lig][pipeIdx[1]].getW()); // OK

					}
				}
				// Debug
	//			cudaDeviceSynchronize();
	//			size_t bufferSize = d_potLig[pipeIdx[1]].bufferSize();
	//			WorkerBuffer<REAL> h_potLig(4,bufferSize);
	//			WorkerBuffer<REAL> h_potRec(4,bufferSize);
	//			size_t cpySize = h_potLig.bufferSize()*sizeof(REAL);
	//			size_t cpySizeRec = h_potRec.bufferSize()*sizeof(REAL);
	//			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[1]].getX(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[1]].getY(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[1]].getW(), cpySize, cudaMemcpyDeviceToHost);
	//
	//			cudaMemcpy(h_potRec.getX(),d_potRec[pipeIdx[1]].getX(), cpySizeRec, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potRec.getY(),d_potRec[pipeIdx[1]].getY(), cpySizeRec, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potRec.getZ(),d_potRec[pipeIdx[1]].getZ(), cpySizeRec, cudaMemcpyDeviceToHost);
	//			cudaMemcpy(h_potRec.getW(),d_potRec[pipeIdx[1]].getW(), cpySizeRec, cudaMemcpyDeviceToHost);
	//
	//			for(size_t i = 0; i < bufferSize; ++i) {
	////			for(size_t i = 0; i < 20; ++i) {
	//
	//		//		std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;
	//				std::cout << h_potRec.getX()[i] << " " << h_potRec.getY()[i] << " " << h_potRec.getZ()[i] << " " << h_potRec.getW()[i] << std::endl;
	//			}
	//			exit(EXIT_SUCCESS);
			}//end loop over ligands
			/* Device: Signal event when force and energy calc. has completed */
			cudaVerify(cudaEventRecord(events[3+pipeIdx[1]], streams[2]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][1] = true;
		}
	}

	void S2_reduce() {
		/* check if stage 1 was executed in last iteration */
		if (predicates[pipeIdx[1]][1] == true) {
			constexpr unsigned stageId = 2;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

			/* Device: Wait for completion of PotForce calc. to complete */
			cudaVerify(cudaStreamWaitEvent(streams[3], events[3+pipeIdx[0]], 0));

			for(unsigned int lig = 0; lig < Common_MB_Modes::numLigands; lig++){
				deviceReduce(
					blockSizeReduce,
					it->size(),
					stageResc.rec->numAtoms,	stageResc.lig[lig]->numAtoms,
					stageResc.rec->numModes, 	stageResc.lig[lig]->numModes,
					d_dof[pipeIdxDof[2]].get(0),
					lig,
					stageResc.rec->xModes, 		stageResc.rec->yModes, 		stageResc.rec->zModes,
					stageResc.lig[lig]->xModes, stageResc.lig[lig]->yModes, stageResc.lig[lig]->zModes,
					stageResc.lig[lig]->xPos,	stageResc.lig[lig]->yPos, 	stageResc.lig[lig]->zPos,
					d_potRec[pipeIdx[0]].getX(),
					d_potRec[pipeIdx[0]].getY(),
					d_potRec[pipeIdx[0]].getZ(),
					d_potLig[lig][pipeIdx[0]].getX(),
					d_potLig[lig][pipeIdx[0]].getY(),
					d_potLig[lig][pipeIdx[0]].getZ(),
					d_potLig[lig][pipeIdx[0]].getW(),
					d_res[pipeIdx[0]].get(0),
					streams[3]);

			cudaDeviceSynchronize();
			unsigned numDofs = it->size();
			WorkerBuffer<REAL> h_potLig(Common_MB_Modes::numLigands,23*numDofs);
			size_t cpySize = h_potLig.bufferSize()*sizeof(REAL);
			cudaMemcpy(h_potLig.get(lig),d_res[pipeIdx[0]].get(lig), cpySize, cudaMemcpyDeviceToHost);
			}
//			for(size_t i = 0; i < numDofs; ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				REAL x = h_potLig.get(0)[i*13 + 0];
//				REAL y = h_potLig.get(0)[i*13 + 1];
//				REAL z = h_potLig.get(0)[i*13 + 2];
//				REAL E = h_potLig.get(0)[i*13 + 3];
//				std::cout << x << " " << y << " " << z << " " << E << std::endl;
//			}
//			for(size_t i = 0; i < 90; ++i) {
//								std::cout <<i <<" "<< h_potLig.get(0)[ i] << std::endl;
//							}
//
//			exit(EXIT_SUCCESS);

			/* Device: Signal event when reduction has completed */
			cudaVerify(cudaEventRecord(events[5+pipeIdx[0]], streams[3]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][2] = true;

		}
	}

	void S3_copyD2H() {
		/* check if stage 2 was executed in last iteration */
		if (predicates[pipeIdx[1]][2])
		{
			constexpr unsigned stageId = 3;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

			/* Device: Wait for completion of reduction calc. to complete */
			cudaVerify(cudaStreamWaitEvent(streams[1], events[5 + pipeIdx[1]], 0));

			/* copy results to host */
			cudaVerify(cudaMemcpyAsync(h_res[pipeIdx[1]].get(0), d_res[pipeIdx[1]].get(0), dofSize*it->size()*sizeof(REAL),
					cudaMemcpyDeviceToHost, streams[1]));

			/* Device: Signal event when transfere has completed */
			cudaVerify(cudaEventRecord(events[1], streams[1]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][3] = true;
		}
	}

	void S4_finalReduce() {
		/* check if stage 3 was executed in last iteration */
		if (predicates[pipeIdx[1]][3] == true)
		{
			constexpr unsigned stageId = 4;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

			/* Host: Wait for completion of data transfer to complete */
			cudaVerify(cudaEventSynchronize(events[1]));

			nvtxRangePushA("Host");
			for(unsigned int lig = 0; lig < Common_MB_Modes::numLigands; lig++){
				h_finalReduce(
					it->size(),
					it->inputBuffer(),
					lig,
					stageResc.rec->modeForce,
					stageResc.lig[lig]->modeForce,
					stageResc.rec->numModes,
					stageResc.lig[lig]->numModes,
					h_res[pipeIdx[0]].get(lig),
					it->resultBuffer());
			}
			nvtxRangePop();

			/* Signal that result is in buffer */
			it->setProcessed();

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][4] = true;

			/* signal that one item has been passed the last stage */
			signalItemPassedLastStage();
		}
	}

	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof[3];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoLig[LIGANDS_MAX_NUMBER];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig[LIGANDS_MAX_NUMBER];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potRec[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig[LIGANDS_MAX_NUMBER ][2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_res[2];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_res[2];

	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	size_t dofSize = 13 + 5 + 5;
	size_t blockSizeReduce = 0;
	static constexpr unsigned numStages = 5;

	bool predicates[2][numStages];
	RingArray<StageResource> stagesMngt;
	unsigned pipeIdx[2];
	unsigned pipeIdxDof[3];

	int numItemsInPipe; /** number of items in the pipeline */

	cudaStream_t streams[4]; /** cuda streams */

	cudaEvent_t events[7];   /** cuda events */



};

template<typename REAL>
auto GPUEnergyService6D_MB_Modes<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		(void)numWorkers;

		std::vector<id_t> ids = {common->gridIdRec, common->recId, common->tableId};
		for(int lig = 0; lig < Common_MB_Modes::numLigands; lig++){
			ids.push_back(common->gridIdLig[lig]);
		}
		for(int lig = 0; lig < Common_MB_Modes::numLigands; lig++){
			ids.push_back(common->ligId[lig]);
		}
		return this->_dataMng->getCommonDeviceIds(ids);
	};
	return fncObj;
}

template<typename REAL>
auto GPUEnergyService6D_MB_Modes<REAL>::createItemProcessor() -> itemProcessor_t {

	std::shared_ptr<Private> p = std::make_shared<Private>();
	deviceId_t deviceId = _workerId++;
	itemProcessor_t fncObj = [this, deviceId, p] (workItem_t* item) -> bool {

		/* Set the device to work with */
		cudaVerify(cudaSetDevice(deviceId));

		/* reset the predicates for the actual iteration*/
		p->resetPrediacatesForIteration();

		p->configureDevice();

		if (item != nullptr) {

			auto dI = createStageResource(item, deviceId);

			const auto itemSize = item->size();
			assert(itemSize > 0);
			const unsigned& numAtomsRec = dI.rec->numAtoms;
			unsigned ligModes = 0;

			p->resizeRecBuffersIfRequired(itemSize, numAtomsRec, Common_MB_Modes::numLigands);
			for( int lig = 0; lig < Common_MB_Modes::numLigands; lig++){
				const unsigned& numAtomsLig = dI.lig[lig]->numAtoms;
				p->addItemAndLigandSize(dI);

				p->resizeLigBuffersIfRequired(itemSize, numAtomsLig, Common_MB_Modes::numLigands, lig );
				ligModes += dI.lig[lig]->numModes;
			}

			const unsigned& dofSize = Common_MB_Modes::numLigands * (13 + ligModes + dI.rec->numModes);
			p->resizeDOFBuffersIfRequired(itemSize, dofSize, Common_MB_Modes::numLigands);
		} else {
			p->stagesMngt.rotate();
		}

		p->S0_copyH2D();
		p->S1_transform_potForce();
		p->S2_reduce();
		p->S3_copyD2H();
		p->S4_finalReduce();

		p->swapBuffers();

		return !(p->pipelineEmpty());

	};

	return fncObj;
}

} // namespace as

#endif

#endif /* SRC_SERVICE_GPUENERGYSERVICE6D_MB_MODES_TPP_ */
