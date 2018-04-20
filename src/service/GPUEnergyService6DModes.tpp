/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */
/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICE6DMODES_TPP_
#define SRC_SERVICE_GPUENERGYSERVICE6DMODES_TPP_

#ifdef CUDA

#include <nvToolsExt.h>
#include "GPUEnergyService6DModes.h"

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

#include "transform_modes.h"
#include "interpolation.h"
#include "neighborlist_modes.h"
#include "reduction_modes.h"

#include "macros.h"



namespace as {

template<typename REAL>
GPUEnergyService6DModes<REAL>::GPUEnergyService6DModes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds) :
	GPUEnergyService<Types_6D_Modes<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds)
{}

template<typename REAL>
struct GPUEnergyService6DModes<REAL>::StageResource {
private:
	using workItem_t = typename GPUEnergyService6DModes<REAL>::workItem_t;
public:
	d_GridUnion<REAL> gridRec;
	d_GridUnion<REAL> gridLig;
	d_Protein<REAL>* rec;
	d_Protein<REAL>* lig;
	d_ParamTable<REAL>* table;
	SimParam<REAL>* simParam;
	workItem_t* item;
};

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createStageResource(workItem_t* item, unsigned const& deviceId) -> StageResource {
	/* item pointers */
//			const auto dofs = item->inputBuffer();
	const auto common = item->common();
//			auto results = item->resultBuffer();

	/* get DataItem pointers */
	auto gridRec = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridIdRec, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
	assert(gridRec != nullptr);

	auto gridLig = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridIdLig, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
	assert(gridLig != nullptr);

	auto rec = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->recId, deviceId)).get();
	assert(rec != nullptr);

	auto lig = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->ligId, deviceId)).get();
	assert(lig != nullptr);

	auto table = std::dynamic_pointer_cast<DeviceParamTable<REAL>>(this->_dataMng->get(common->tableId, deviceId)).get();
	assert(table != nullptr);

	auto simParam = std::dynamic_pointer_cast<SimParam<REAL>>(this->_dataMng->get(common->paramsId)).get();
	assert(simParam != nullptr);

	StageResource stageResource;
	stageResource.gridRec 	= gridRec->getDesc();
	stageResource.gridLig 	= gridLig->getDesc();
	stageResource.lig 		= &lig->desc;
	stageResource.rec 		= &rec->desc;
	stageResource.table 	= &table->desc;
	stageResource.simParam = simParam;
	stageResource.item 	= item;

	return stageResource;
}

template<typename REAL>
class GPUEnergyService6DModes<REAL>::Private {
	using dof_t = typename GPUEnergyService6DModes<REAL>::input_t;
	using workItem_t = typename GPUEnergyService6DModes<REAL>::workItem_t;

public:

	Private() : stagesMngt(numStages), pipeIdx{0,1}, pipeIdxDof{0,1,2},numItemsInPipe(0) {
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
	void allocateBuffer( size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numAtomsLig, size_t const& dofSizeRec, size_t const& dofSizeLig ) {
		const size_t atomBufferSizeRec = numDOFs*numAtomsRec;
		const size_t atomBufferSizeLig = numDOFs*numAtomsLig;
		for (int i = 0; i < 2; ++i) {
			d_potRec[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSizeRec));
			d_potLig[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSizeLig));
			d_resRec[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSizeRec*numDOFs));
			d_resLig[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSizeLig*numDOFs));
			h_resRec[i] = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSizeRec*numDOFs));
			h_resLig[i] = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSizeLig*numDOFs));
		}
		for (int i = 0 ; i < 4; ++i){
			d_dof[i]    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
		}
		d_defoRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeRec));
		d_defoLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeLig));
		d_trafoRec = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeRec));
		d_trafoLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeLig));
	}

	size_t bufferSize() const {
		return d_trafoLig.bufferSize();
	}

	void addItemAndLigandSize(StageResource const& resc) {
		++numItemsInPipe;
		predicates[pipeIdx[0]][0] = true;
		stagesMngt.push(resc);
	}

	void resizeBuffersIfRequired(size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numAtomsLig, size_t const& dofSizeRec, size_t const& dofSizeLig ) {
		if (numDOFs*numAtomsLig > d_trafoLig.bufferSize() || numDOFs*numAtomsRec > d_trafoRec.bufferSize()) {
			allocateBuffer(numDOFs, numAtomsRec, numAtomsLig, dofSizeRec, dofSizeLig );
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


			cudaVerify(cudaStreamWaitEvent(streams[0], events[2], 0));
			cudaVerify(cudaMemcpyAsync(d_dof[pipeIdxDof[0]].get(0), it->inputBuffer(),
					it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, streams[0]));
			cudaVerify(cudaEventRecord(events[0], streams[0]));

		}
	}

	size_t getSharedMemSize(int const& id){
		cudaDeviceProp deviceProp;
		cudaVerify(cudaGetDeviceProperties(&deviceProp, id));
		return deviceProp.sharedMemPerBlock;
	}

	size_t getMaxBockSize( int const& id, unsigned const dofSize ){
		size_t sharedMem = getSharedMemSize( id );
		size_t pow2 = 2;
		while (pow2* dofSize *sizeof(REAL) < sharedMem) {
			pow2 *= 2;
		}
		size_t blockSize = pow2 / 2;
		return blockSize;
	}

	void configureDevice() {

		if ( blockSizeReduceRec == 0 ) {
			int id;
			cudaVerify(cudaGetDevice(&id));
			blockSizeReduceRec = getMaxBockSize( id, dofSizeRec );
		}
		if ( blockSizeReduceLig == 0 ) {
			int id;
			cudaVerify(cudaGetDevice(&id));
			blockSizeReduceLig = getMaxBockSize( id, dofSizeLig );
		}

	}

	void S1_transform_potForce() {
		/* check if stage 0 was executed in last iteration */
		if (predicates[pipeIdx[1]][0])
		{
			constexpr unsigned stageId = 1;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

			const unsigned numElRec = it->size()*stageResc.rec->numAtoms;
			const unsigned numElLig = it->size()*stageResc.lig->numAtoms;

			assert( numElRec <= d_trafoRec.bufferSize() );
			assert( numElLig <= d_trafoLig.bufferSize() );

			/* Perform cuda kernel calls */
			size_t gridSizeRec = ( numElRec + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
			size_t gridSizeLig = ( numElLig + BLSZ_TRAFO - 1) / BLSZ_TRAFO;

			/* Device: Wait for completion of copyH2D of DOFs to complete */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[0], 0));
			stageResc.rec->numModes = 5;
			stageResc.lig->numModes = 5;

			d_DOFPos_rec<REAL,  true>(
				BLSZ_TRAFO,
				gridSizeRec,
				streams[2],
				stageResc.rec,
				d_dof[pipeIdxDof[1]].get(0),
				it->size(),
				d_defoRec.getX(),
				d_defoRec.getY(),
				d_defoRec.getZ(),
				d_trafoRec.getX(),
				d_trafoRec.getY(),
				d_trafoRec.getZ()
				);

			d_DOFPos_lig<REAL,true>(
				BLSZ_TRAFO,
				gridSizeLig,
				streams[2],
				stageResc.lig,
				d_dof[pipeIdxDof[1]].get(0),
				it->size(),
				d_defoLig.getX(),
				d_defoLig.getY(),
				d_defoLig.getZ(),
				d_trafoLig.getX(),
				d_trafoLig.getY(),
				d_trafoLig.getZ()
				);

////			 DEBUG
//			cudaDeviceSynchronize();
//			size_t bufferSize = d_defoLig.bufferSize();
//			WorkerBuffer<REAL> h_defoLig(3,bufferSize);
//			size_t cpySize = h_defoLig.bufferSize()*sizeof(REAL);
//
//			std::cout << "bufferSize: " << bufferSize << " cpySize: " << cpySize << std::endl;
//			cudaMemcpy(h_defoLig.getX(),d_defoLig.getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_defoLig.getY(),d_defoLig.getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_defoLig.getZ(),d_defoLig.getZ(), cpySize, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < bufferSize; ++i) {
//				std::cout << h_defoLig.getX()[i] << " " << h_defoLig.getY()[i] << " " << h_defoLig.getZ()[i] << std::endl;
//			}
			//exit(EXIT_SUCCESS);

			/* Device: Signal event when transformation has completed */
			cudaVerify(cudaEventRecord(events[2], streams[2]));
			/* Device: Wait for completion of reduction of the previous round */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[5+pipeIdx[1]], 0));

			/* Perform cuda kernel calls */
			gridSizeRec = ( numElRec + BLSZ_INTRPL - 1) / BLSZ_INTRPL;
			gridSizeLig = ( numElLig + BLSZ_INTRPL - 1) / BLSZ_INTRPL;

			d_potForce (
				BLSZ_INTRPL,
				gridSizeLig,
				streams[2],
				stageResc.gridRec.inner,
				stageResc.gridRec.outer,
				*stageResc.lig,
				it->size(),
				d_trafoLig.getX(),
				d_trafoLig.getY(),
				d_trafoLig.getZ(),
				d_potLig[pipeIdx[1]].getX(),
				d_potLig[pipeIdx[1]].getY(),
				d_potLig[pipeIdx[1]].getZ(),
				d_potLig[pipeIdx[1]].getW()); // OK

			d_potForce (
				BLSZ_INTRPL,
				gridSizeRec,
				streams[2],
				stageResc.gridLig.inner,
				stageResc.gridLig.outer,
				*stageResc.rec,
				it->size(),
				d_trafoRec.getX(),
				d_trafoRec.getY(),
				d_trafoRec.getZ(),
				d_potRec[pipeIdx[1]].getX(),
				d_potRec[pipeIdx[1]].getY(),
				d_potRec[pipeIdx[1]].getZ(),
				d_potRec[pipeIdx[1]].getW()); // OK



			//std::cout << sumx <<  " " << sumy << " " << sumz;
			//exit(EXIT_SUCCESS);

			// Debug


			d_NLPotForce(
				BLSZ_INTRPL,
				gridSizeLig,
				streams[2],
				stageResc.gridRec.NL,
				*stageResc.rec,
				*stageResc.lig,
				*stageResc.table,
				*stageResc.simParam,
				it->size(),
				d_defoRec.getX(),
				d_defoRec.getY(),
				d_defoRec.getZ(),
				d_trafoLig.getX(),
				d_trafoLig.getY(),
				d_trafoLig.getZ(),
				d_potLig[pipeIdx[1]].getX(),
				d_potLig[pipeIdx[1]].getY(),
				d_potLig[pipeIdx[1]].getZ(),
				d_potLig[pipeIdx[1]].getW()); // OK


//
//			//get nl forces on receptor
			d_NLPotForce(
				BLSZ_INTRPL,
				gridSizeRec,
				streams[2],
				stageResc.gridLig.NL,
				*stageResc.lig,
				*stageResc.rec,
				*stageResc.table,
				*stageResc.simParam,
				it->size(),
				d_defoLig.getX(),
				d_defoLig.getY(),
				d_defoLig.getZ(),
				d_trafoRec.getX(),
				d_trafoRec.getY(),
				d_trafoRec.getZ(),
				d_potRec[pipeIdx[1]].getX(),
				d_potRec[pipeIdx[1]].getY(),
				d_potRec[pipeIdx[1]].getZ(),
				d_potRec[pipeIdx[1]].getW()
				); // OK
//
//
//			cudaDeviceSynchronize();

//			WorkerBuffer<REAL> h_potLig(4,3*stageResc.lig->numAtoms);
//			size_t cpySize = 3*stageResc.lig->numAtoms*sizeof(REAL);
//			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[1]].getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[1]].getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[1]].getW(), cpySize, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < 3*stageResc.lig->numAtoms; ++i) {
//				//			for(size_t i = 0; i < 20; ++i) {
//				std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;

//			}


			d_rotateForces(
				BLSZ_INTRPL,
				gridSizeRec,
				streams[2],
				d_potRec[pipeIdx[1]].getX(),
				d_potRec[pipeIdx[1]].getY(),
				d_potRec[pipeIdx[1]].getZ(),
				d_dof[pipeIdxDof[1]].get(0),
				stageResc.rec->numAtoms,
				it->size()
				);
			//IP.d_NLPotForce<false>(it->devLocGridId(), it->devLocRecId(), it->devLocRecId(),it->size(),
 		//	&d_trafoRec, d_potRec[pipeIdx[1]],streams[2]);


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


//			cudaDeviceSynchronize();
//			size_t bufferSize = d_potLig[pipeIdx[0]].bufferSize();
//			WorkerBuffer<REAL> h_potLig(4,bufferSize);
//			size_t cpySize = h_potLig.bufferSize()*sizeof(REAL);
//			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[0]].getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[0]].getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[0]].getZ(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[0]].getW(), cpySize, cudaMemcpyDeviceToHost);
//			float esum = 0 ;
//			for(size_t i = 0; i < 875; ++i) {
////			for(size_t i = 0; i < 20;  ++i) {
//				esum += h_potLig.getW()[i];
//				std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;
//			}
//			std::cout << esum ;
//			exit(EXIT_SUCCESS);

			deviceReduce<REAL, 0, true>(
				blockSizeReduceRec,
				it->size(),
				stageResc.rec,
				d_dof[pipeIdxDof[2]].get(0),
				stageResc.rec->xPos, stageResc.rec->yPos, stageResc.rec->zPos,
				d_potRec[pipeIdx[0]].getX(), d_potRec[pipeIdx[0]].getY(), d_potRec[pipeIdx[0]].getZ(),
				d_potRec[pipeIdx[0]].getW(),
				d_resRec[pipeIdx[0]].get(0),
				streams[3]);

			deviceReduce<REAL, 1, true>(
				blockSizeReduceLig,
				it->size(),
				stageResc.lig,
				d_dof[pipeIdxDof[2]].get(0),
				stageResc.lig->xPos, stageResc.lig->yPos, stageResc.lig->zPos,
				d_potLig[pipeIdx[0]].getX(), d_potLig[pipeIdx[0]].getY(), d_potLig[pipeIdx[0]].getZ(),
				d_potLig[pipeIdx[0]].getW(),
				d_resLig[pipeIdx[0]].get(0),
				streams[3]);



//			cudaDeviceSynchronize();
//			unsigned numDofs = it->size();
//
//			WorkerBuffer<REAL> h_potLig(1,(dofSizeLig)*numDofs);
//			size_t cpySize = h_potLig.bufferSize()*sizeof(REAL);
//			cudaMemcpy(h_potLig.get(0),d_resLig[pipeIdx[0]].get(0), cpySize, cudaMemcpyDeviceToHost);
//
/////			for(size_t i = 0; i < numDofs; ++i) {
//			for(size_t i = 0; i < 1; ++i) {
//				REAL x = h_potLig.get(0)[i*dofSizeLig + 0];
//				REAL y = h_potLig.get(0)[i*dofSizeLig + 1];
//				REAL z = h_potLig.get(0)[i*dofSizeLig + 2];
//				REAL E = h_potLig.get(0)[i*dofSizeLig + 3];
//				REAL lm1 = h_potLig.get(0)[i*dofSizeLig + 13];
//				REAL lm2 = h_potLig.get(0)[i*dofSizeLig + 14];
//				REAL lm3 = h_potLig.get(0)[i*dofSizeLig + 15];
//				REAL lm4 = h_potLig.get(0)[i*dofSizeLig + 16];
//				REAL lm5 = h_potLig.get(0)[i*dofSizeLig + 17];
//
//				std::cout << x << " " << y << " " << z << " " << E <<" " << lm1 <<" " << lm2<<" " << lm3 <<" " << lm4<<" " << lm5 << std::endl;
//			}
			//exit(EXIT_SUCCESS);

			/* Device: Signal event when reduction has completed */
			cudaVerify(cudaEventRecord(events[5+pipeIdx[0]], streams[3]));
			//std::cout << " endhello check " << std::endl;
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
			cudaVerify(cudaMemcpyAsync(h_resRec[pipeIdx[1]].get(0), d_resRec[pipeIdx[1]].get(0), dofSizeRec*it->size()*sizeof(REAL),
					cudaMemcpyDeviceToHost, streams[1]));
			cudaVerify(cudaMemcpyAsync(h_resLig[pipeIdx[1]].get(0), d_resLig[pipeIdx[1]].get(0), dofSizeLig*it->size()*sizeof(REAL),
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
//			std::cout << std::endl;
//
//			for ( int i= 0; i < 3; i++){
//				std::cout << " new calculation" << std::endl;
//				std::cout << " " << h_resLig[pipeIdx[0]].bufferSize()<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i+1]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i + 2]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i + 3]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i + 13]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i + 14]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i + 15]<< std::endl;
//							std::cout << " " << h_resLig[pipeIdx[0]].get(0)[18*i + 16]<< std::endl;
//						}
////			for ( int i= 13; i < 18; i++){
////				std::cout << " " << h_resLig[pipeIdx[1]].get(0)[i];
////			}
//			std::cout << std::endl;

			nvtxRangePushA("Host");

			h_finalReduce< REAL, 1, true>(
				it->size(),
				stageResc.lig,
				it->inputBuffer(),
				h_resLig[pipeIdx[0]].get(0),
				it->resultBuffer());
				nvtxRangePop();

			h_finalReduce< REAL, 0, true>(
				it->size(),
				stageResc.rec,
				it->inputBuffer(),
				h_resRec[pipeIdx[0]].get(0),
				it->resultBuffer());
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
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoRec;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potRec[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_resRec[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_resLig[2];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_resRec[2];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_resLig[2];


	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	size_t blockSizeReduceRec = 0;
	size_t blockSizeReduceLig = 0;
	size_t dofSizeRec = 13 + Common_Modes::numModesRec;
	size_t dofSizeLig = 13 + Common_Modes::numModesLig;
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
auto GPUEnergyService6DModes<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		(void)numWorkers;
		std::vector<id_t> ids = {common->gridIdRec, common->gridIdLig, common->ligId, common->recId, common->tableId};
		return this->_dataMng->getCommonDeviceIds(ids);
	};
	return fncObj;
}

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createItemProcessor() -> itemProcessor_t {

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
			const unsigned& numAtomsLig = dI.lig->numAtoms;
			const unsigned& dofSizeRec = 13 + Common_Modes::numModesRec;
			const unsigned& dofSizeLig= 13 + Common_Modes::numModesLig;

			p->addItemAndLigandSize(dI);
			p->resizeBuffersIfRequired(itemSize, numAtomsRec, numAtomsLig, dofSizeRec, dofSizeLig );

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

#endif /* SRC_SERVICE_GPUENERGYSERVICE6D_TPP_ */
