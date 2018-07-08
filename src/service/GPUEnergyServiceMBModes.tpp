/*
 * GPU_MB_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */
/*
 * GPU_MB_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICEMBMODES_TPP_
#define SRC_SERVICE_GPUENERGYSERVICEMBMODES_TPP_

#ifdef CUDA

#include <nvToolsExt.h>
#include "GPUEnergyServiceMBModes.h"

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
#include "scoring_kernel.h"
#include "macros.h"
#include <iostream>
#include "ThreadSafeQueue.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <mutex>
#include "Types_MB_Modes.h"
namespace as {

template<typename REAL>
GPUEnergyServiceMBModes<REAL>::GPUEnergyServiceMBModes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds) :
	GPUEnergyService<Types_MB_Modes<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds)
{}

template<typename REAL>
struct GPUEnergyServiceMBModes<REAL>::StageResource {
private:
	using workItem_t = typename GPUEnergyServiceMBModes<REAL>::workItem_t;
public:
	d_GridUnion<REAL> grids[NUM_MAX_PROTEIN];
	d_Protein<REAL>* proteins[NUM_MAX_PROTEIN];
	d_ParamTable<REAL>* table;
	SimParam<REAL>* simParam;
	workItem_t* item;
};

template<typename REAL>
auto GPUEnergyServiceMBModes<REAL>::createStageResource(workItem_t* item, unsigned const& deviceId) -> StageResource {
	/* item pointers */
//			const auto dofs = item->inputBuffer();
	const auto common = item->common();
	StageResource stageResource;
//			auto results = item->resultBuffer();

	for( int id = 0; id < common->numProteins; ++id){
		auto protConfig = common->proteins[id];
		auto grid = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(protConfig.gridId, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(grid != nullptr);
		auto protein = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(protConfig.proteinId, deviceId)).get();
		assert(protein != nullptr);

		stageResource.grids[id] 		= grid->getDesc();
		stageResource.proteins[id] 		= &protein->desc;

	}
	/* get DataItem pointers */


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
class GPUEnergyServiceMBModes<REAL>::Private {
	using dof_t = typename GPUEnergyServiceMBModes<REAL>::input_t;
	using workItem_t = typename GPUEnergyServiceMBModes<REAL>::workItem_t;

public:

	Private() : numItemsInPipe(0) {
		for ( unsigned id_protein = 0; id_protein < Common_MB_Modes::numProteins; ++id_protein){
			dofSize[id_protein] = 13 + Common_MB_Modes::numModes[id_protein];
		}

			CUDA_CHECK(cudaStreamCreate(&stream));
	}
	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBuffer( size_t const& numDOFs, size_t const& numAtoms ,size_t const& dofSize, unsigned const id_protein ) {
		const size_t atomBufferSize = numDOFs*numAtoms;

		d_defo[id_protein] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSize));

		d_pot[0][id_protein] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
		d_pot[1][id_protein] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
		d_res[id_protein] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSize*numDOFs));
		h_res[id_protein] = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSize*numDOFs));
		d_dof    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
		h_dofSupport = std::move(WorkerBuffer<dofMB<REAL>,HostPinnedAllocator<dofMB<REAL>>>(1,numDOFs*Common_MB_Modes::numProteins*Common_MB_Modes::numProteins));
		d_dofSupport = std::move(WorkerBuffer<dofMB<REAL>,DeviceAllocator<dofMB<REAL>>>(1,numDOFs*Common_MB_Modes::numProteins*Common_MB_Modes::numProteins));
		for ( int i = 0; i< Common_MB_Modes::numProteins; ++i)
		{
			d_trafo[id_protein][i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSize));
		}

	}

	size_t bufferSize( unsigned const id_protein ) const {
		return d_trafo[id_protein][0].bufferSize(  );
	}

	void addItemAndLigandSize(StageResource const& resc) {
		_lock.lock();
		_resource = resc;
		++numItemsInPipe;
		_lock.unlock();
	}

	void resizeBuffersIfRequired(size_t const& numDOFs, size_t const& numAtoms, size_t const& dofSize,  unsigned const id_protein) {
		if (numDOFs*numAtoms > d_trafo[id_protein][0].bufferSize()) {
			allocateBuffer(numDOFs, numAtoms,  dofSize, id_protein );
		}
	}

	bool pipelineEmpty() const {
		return numItemsInPipe == 0;
	}

	void signalItemPassedLastStage() {
		_lock.lock();
		--numItemsInPipe;
		_lock.unlock();
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
		for( unsigned id_protein = 0; id_protein < Common_MB_Modes::numProteins; ++id_protein)
		{
			if ( blockSizeReduceRec == 0 ) {
				int id;
				cudaVerify(cudaGetDevice(&id));
				blockSizeReduceRec = getMaxBockSize( id, dofSize[id_protein]);
			}
		}
	}

	void score( ) {

		auto const& stageResc = _resource;
		auto* const it = stageResc.item;
		const auto common = it->common();
		std::vector<Vec3<REAL>> pivots;
		for( int i = 0; i< common->numProteins; ++i){

			if ( !common->proteins[i].centered ){
				//std::cout<< "sdgasgasg" << common->proteins[i].pivot << std::endl;
				pivots.push_back(Vec3<REAL>(
						common->proteins[i].pivot.x,
						common->proteins[i].pivot.y,
						common->proteins[i].pivot.z
						));

			}
			else{
				pivots.push_back(Vec3<REAL>(0.0f));
			}
		}
		cudaVerify(cudaMemcpyAsync(d_dof.get(0), it->inputBuffer(),
						it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, stream));
		getNewDof<REAL>( it->size(), common->numProteins,  it->inputBuffer(),pivots, h_dofSupport.get(0));
		cudaVerify(cudaMemcpyAsync(d_dofSupport.get(0), h_dofSupport.get(0),
						it->size()*common->numProteins*common->numProteins*sizeof(dofMB<REAL>), cudaMemcpyHostToDevice, stream));

		for( int id_centerProtein = 0; id_centerProtein < common->numProteins ; ++id_centerProtein){
			cudaMemsetAsync(d_pot[0][id_centerProtein].getX(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[0][id_centerProtein].getY(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[0][id_centerProtein].getZ(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[0][id_centerProtein].getW(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[1][id_centerProtein].getX(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[1][id_centerProtein].getY(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[1][id_centerProtein].getZ(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
			cudaMemsetAsync(d_pot[1][id_centerProtein].getW(), 0, it->size()*stageResc.proteins[id_centerProtein]->numAtoms* sizeof(REAL), stream);
		}

		for ( unsigned id_centerProtein = 0; id_centerProtein< common->numProteins; ++id_centerProtein)
		{
			unsigned const numAtomsCenter = it->size()*stageResc.proteins[id_centerProtein]->numAtoms;
			assert( numAtomsCenter <= d_trafo[id_centerProtein][0].bufferSize() );
			size_t gridSizeCenter = ( numAtomsCenter + BLSZ_TRAFO - 1) / BLSZ_TRAFO;

			//deform needs unpivotized and transform needs pivotized coordinates
			 d_deform(
					 BLSZ_TRAFO,
					gridSizeCenter,
					stream,
					*stageResc.proteins[id_centerProtein],
					d_dof.get(0),
					it->size(),
					id_centerProtein,
					stageResc.proteins[id_centerProtein]->xPos,
					stageResc.proteins[id_centerProtein]->yPos,
					stageResc.proteins[id_centerProtein]->zPos,
					d_defo[id_centerProtein].getX(),
					d_defo[id_centerProtein].getY(),
					d_defo[id_centerProtein].getZ()
					);

		}
		for ( unsigned id_centerProtein = 0; id_centerProtein< common->numProteins; ++id_centerProtein)
				{
					unsigned const numAtomsCenter = it->size()*stageResc.proteins[id_centerProtein]->numAtoms;
					assert( numAtomsCenter <= d_trafo[id_centerProtein][0].bufferSize() );
					size_t gridSizeCenter = ( numAtomsCenter + BLSZ_TRAFO - 1) / BLSZ_TRAFO;

			 for ( unsigned id_partnerProtein = 0; id_partnerProtein< common->numProteins; ++id_partnerProtein)
			{
				 if(id_partnerProtein != id_centerProtein){

					 d_transform_new(
						BLSZ_TRAFO,
						gridSizeCenter,
						stream,
						d_dofSupport.get(0),
						it->size(),
						stageResc.proteins[id_centerProtein]->numAtoms,
						common->numProteins,
						id_centerProtein,
						id_partnerProtein,
						d_defo[id_centerProtein].getX(),
						d_defo[id_centerProtein].getY(),
						d_defo[id_centerProtein].getZ(),
						d_trafo[id_centerProtein][id_partnerProtein].getX(),
						d_trafo[id_centerProtein][id_partnerProtein].getY(),
						d_trafo[id_centerProtein][id_partnerProtein].getZ()
							);				 }
			 }

		}


		for ( unsigned id_centerProtein = 0; id_centerProtein< common->numProteins; ++id_centerProtein)
		{
			unsigned const numAtomsCenter = it->size()*stageResc.proteins[id_centerProtein]->numAtoms;
			assert( numAtomsCenter <= d_trafo[id_centerProtein][0].bufferSize() );
			size_t gridSizeCenter = ( numAtomsCenter + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
			for ( unsigned id_partnerProtein = 0; id_partnerProtein< common->numProteins; ++id_partnerProtein)
			{
				if (common->radius_cutoff > 10000)
				{
					if(id_partnerProtein != id_centerProtein){
					d_potForce(
						BLSZ_TRAFO,
						gridSizeCenter,
						stream,
						stageResc.grids[id_partnerProtein].inner,
						stageResc.grids[id_partnerProtein].outer,
						*stageResc.proteins[id_centerProtein],
						it->size(),
						d_trafo[id_centerProtein][id_partnerProtein].getX(),
						d_trafo[id_centerProtein][id_partnerProtein].getY(),
						d_trafo[id_centerProtein][id_partnerProtein].getZ(),
						d_pot[0][id_centerProtein].getX(),
						d_pot[0][id_centerProtein].getY(),
						d_pot[0][id_centerProtein].getZ(),
						d_pot[0][id_centerProtein].getW());


				d_NLPotForce<REAL, dof_t>(
						BLSZ_INTRPL,
						gridSizeCenter,
						stream,
						stageResc.grids[id_partnerProtein].NL,
						*stageResc.proteins[id_partnerProtein],
						*stageResc.proteins[id_centerProtein],
						*stageResc.table,
						common->radius_cutoff,
						*stageResc.simParam,
						it->size(),
						d_defo[id_partnerProtein].getX(),
						d_defo[id_partnerProtein].getY(),
						d_defo[id_partnerProtein].getZ(),
						d_trafo[id_centerProtein][id_partnerProtein].getX(),
						d_trafo[id_centerProtein][id_partnerProtein].getY(),
						d_trafo[id_centerProtein][id_partnerProtein].getZ(),
						d_pot[0][id_centerProtein].getX(),
						d_pot[0][id_centerProtein].getY(),
						d_pot[0][id_centerProtein].getZ(),
						d_pot[0][id_centerProtein].getW()
						);
				d_rotateForces(
						BLSZ_INTRPL,
						gridSizeCenter,
						stream,
						id_partnerProtein,
						d_pot[0][id_centerProtein].getX(),
						d_pot[0][id_centerProtein].getY(),
						d_pot[0][id_centerProtein].getZ(),
						d_pot[0][id_centerProtein].getW(),
						d_pot[1][id_centerProtein].getX(),
						d_pot[1][id_centerProtein].getY(),
						d_pot[1][id_centerProtein].getZ(),
						d_pot[1][id_centerProtein].getW(),
						d_dof.get(0),
						stageResc.proteins[id_centerProtein]->numAtoms,
						it->size()
						);

//				cudaDeviceSynchronize();
//				WorkerBuffer<REAL> h_potLig(4,stageResc.proteins[id_centerProtein]->numAtoms);
//				size_t cpySize = stageResc.proteins[id_centerProtein]->numAtoms*sizeof(REAL);
//				//std::cout <<"fx fy fz"<<std::endl;
//				std::cout << "nl Force " << std::endl;
//				cudaMemcpy(h_potLig.getX(),d_pot[id_centerProtein].getX(), cpySize, cudaMemcpyDeviceToHost);
//				cudaMemcpy(h_potLig.getY(),d_pot[id_centerProtein].getY(), cpySize, cudaMemcpyDeviceToHost);
//				cudaMemcpy(h_potLig.getZ(),d_pot[id_centerProtein].getZ(), cpySize, cudaMemcpyDeviceToHost);
//				cudaMemcpy(h_potLig.getW(),d_pot[id_centerProtein].getW(), cpySize, cudaMemcpyDeviceToHost);
//				for(size_t i = 0; i < stageResc.proteins[id_centerProtein]->numAtoms; ++i) {
//					//			for(size_t i = 0; i < 20; ++i) {
//					std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i]<< std::endl;// << " " << h_potLig.getW()[i] ;
//
//				}
//				exit(EXIT_SUCCESS);
					}
				}
			}}

//	cudaDeviceSynchronize();
//	std::cout <<"defoRec g"<<std::endl;
//	size_t bufferSizeDefoRec1 = d_defo[id_partnerProtein].bufferSize();
//	WorkerBuffer<REAL> h_DefoRec(3,bufferSizeDefoRec1);
//	size_t cpySizeDefoRec1 = h_DefoRec.bufferSize()*sizeof(REAL);
//
//	//std::cout << "bufferSize: " << bufferSizeDefoRec << " cpySize: " << cpySizeDefoRec << std::endl;
//	cudaMemcpy(h_DefoRec.getX(),d_defo[id_partnerProtein].getX(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_DefoRec.getY(),d_defo[id_partnerProtein].getY(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_DefoRec.getZ(),d_defo[id_partnerProtein].getZ(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
//	for(size_t i = 0; i < stageResc.proteins[0]->numAtoms; ++i) {
//	std::cout  << std::setprecision(10)<< h_DefoRec.getX()[i] << " " << h_DefoRec.getY()[i] << " " << h_DefoRec.getZ()[i] << std::endl;
//	}
//	std::cout <<"trafoRec g"<<std::endl;
//	size_t bufferSizetrafoRec1 = d_trafo[id_centerProtein][id_partnerProtein].bufferSize();
//	WorkerBuffer<REAL> h_trafoRec(3,bufferSizetrafoRec1);
//	size_t cpySizetrafoRec1 = h_trafoRec.bufferSize()*sizeof(REAL);
//
//	cudaMemcpy(h_trafoRec.getX(),d_trafo[id_centerProtein][id_partnerProtein].getX(), cpySizetrafoRec1, cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_trafoRec.getY(),d_trafo[id_centerProtein][id_partnerProtein].getY(), cpySizetrafoRec1, cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_trafoRec.getZ(),d_trafo[id_centerProtein][id_partnerProtein].getZ(), cpySizetrafoRec1, cudaMemcpyDeviceToHost);
//	for(size_t i = 0; i < stageResc.proteins[1]->numAtoms; ++i) {
//	std::cout  << std::setprecision(10)<< h_trafoRec.getX()[i] << " " << h_trafoRec.getY()[i] << " " << h_trafoRec.getZ()[i] << std::endl;
//	}

	for ( unsigned id_centerProtein = 0; id_centerProtein< common->numProteins; ++id_centerProtein)
				{
			unsigned const numAtomsCenter = it->size()*stageResc.proteins[id_centerProtein]->numAtoms;//
			assert( numAtomsCenter <= d_trafo[id_centerProtein][0].bufferSize() );
			size_t gridSizeCenter = ( numAtomsCenter + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
			deviceReduce<REAL, DOF_MB_Modes<REAL>,1>(
			blockSizeReduceRec,
			it->size(),
			id_centerProtein,
			stageResc.proteins[id_centerProtein],
			d_dof.get(0),
			d_defo[id_centerProtein].getX(), d_defo[id_centerProtein].getY(), d_defo[id_centerProtein].getZ(),
			d_pot[1][id_centerProtein].getX(), d_pot[1][id_centerProtein].getY(), d_pot[1][id_centerProtein].getZ(),
			d_pot[1][id_centerProtein].getW(),
			d_res[id_centerProtein].get(0),
			stream);


			cudaVerify(cudaMemcpyAsync( h_res[id_centerProtein].get(0), d_res[id_centerProtein].get(0), dofSize[id_centerProtein]*it->size()*sizeof(REAL),
					cudaMemcpyDeviceToHost, stream));
			cudaVerify(cudaStreamSynchronize(stream));
			nvtxRangePushA("Host");

			h_finalReduce< REAL>(
				it->size(),
				id_centerProtein,
				stageResc.proteins[id_centerProtein],
				it->inputBuffer(),
				h_res[id_centerProtein].get(0),
				it->resultBuffer());

			nvtxRangePop();
		}
		for ( int i = 0; i < it->size(); ++i){
			it->resultBuffer()[i].E = it->resultBuffer()[i].E / (float)common->numProteins;
		}

//						cudaDeviceSynchronize();
//
//			std::cout <<"defoRec g"<<std::endl;
//			size_t bufferSizeDefoRec1 = d_defoRec[id_stream].bufferSize();
//			WorkerBuffer<REAL> h_DefoRec(3,bufferSizeDefoRec1);
//			size_t cpySizeDefoRec1 = h_DefoRec.bufferSize()*sizeof(REAL);
//
//			//std::cout << "bufferSize: " << bufferSizeDefoRec << " cpySize: " << cpySizeDefoRec << std::endl;
//			cudaMemcpy(h_DefoRec.getX(),d_defoRec[id_stream].getX(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_DefoRec.getY(),d_defoRec[id_stream].getY(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_DefoRec.getZ(),d_defoRec[id_stream].getZ(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < stageResc.rec->numAtoms; ++i) {
//				std::cout  << std::setprecision(10)<< h_DefoRec.getX()[i] << " " << h_DefoRec.getY()[i] << " " << h_DefoRec.getZ()[i] << std::endl;
//			}
////
//
//
//
//
//			std::cout <<"trafo lig "<<std::endl;
//			size_t bufferSizeTrafoLig = d_trafoLig[id_stream].bufferSize();
//
//			WorkerBuffer<REAL> h_trafoLig(3,bufferSizeTrafoLig);
//			size_t cpySizetrafoLig = h_trafoLig.bufferSize()*sizeof(REAL);
////
////			std::cout << "bufferSize: " << bufferSizeDefoLig << " cpySize: " << cpySizeDefoLig << std::endl;
//			cudaMemcpy(h_trafoLig.getX(),d_trafoLig[id_stream].getX(), cpySizetrafoLig, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_trafoLig.getY(),d_trafoLig[id_stream].getY(), cpySizetrafoLig, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_trafoLig.getZ(),d_trafoLig[id_stream].getZ(), cpySizetrafoLig, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
//				std::cout  << std::setprecision(10)<<h_trafoLig.getX()[i] << " " << h_trafoLig.getY()[i] << " " << h_trafoLig.getZ()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);
//
			/* Perform cuda kernel calls */
//			gridSizeRec = ( numElRec + BLSZ_INTRPL - 1) / BLSZ_INTRPL;
//			gridSizeLig = ( numElLig + BLSZ_INTRPL - 1) / BLSZ_INTRPL;

			//std::cout <<"before nl gpu" <<std::endl;
//			cudaDeviceSynchronize();
//			WorkerBuffer<REAL> h_potLig(4,stageResc.lig->numAtoms);
//			size_t cpySize = stageResc.lig->numAtoms*sizeof(REAL);
//			//std::cout <<"fx fy fz"<<std::endl;
//			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[1]].getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[1]].getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[1]].getW(), cpySize, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
//				//			for(size_t i = 0; i < 20; ++i) {
//				std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i]<< std::endl;// << " " << h_potLig.getW()[i] ;
//
//			}






//			std::cout <<"after nl"<< std::endl;
//			cudaDeviceSynchronize();
//			WorkerBuffer<REAL> h_potLig1(4,stageResc.lig->numAtoms);
//			size_t cpySize1 = stageResc.lig->numAtoms*sizeof(REAL);
//			std::cout <<"fx fy fz"<<std::endl;
//			cudaMemcpy(h_potLig1.getX(),d_potLig[pipeIdx[1]].getX(), cpySize1, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig1.getY(),d_potLig[pipeIdx[1]].getY(), cpySize1, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig1.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize1, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig1.getW(),d_potLig[pipeIdx[1]].getW(), cpySize1, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
//				//			for(size_t i = 0; i < 20; ++i) {
//				std::cout << h_potLig1.getX()[i] << " " << h_potLig1.getY()[i] << " " << h_potLig1.getZ()[i]<< std::endl;// << " " << h_potLig1.getW()[i] ;
//
//			}
//			exit(EXIT_SUCCESS);
			//IP.d_NLPotForce<false>(it->devLocGridId(), it->devLocRecId(), it->devLocRecId(),it->size(),
 		//	&d_trafoRec, d_potRec[pipeIdx[1]],streams[2]);



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




			/* copy results to host */



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
//			for ( int i= 13; i < 18; i++){
//				std::cout << " " << h_resLig[pipeIdx[1]].get(0)[i];
//			}
			//std::cout << std::endl;



			/* Signal that result is in buffer */
			it->setProcessed();

			/* signal that this stage was executed within the current iteration */


			/* signal that one item has been passed the last stage */
			signalItemPassedLastStage();
	}
	static unsigned constexpr num_proteins = NUM_MAX_PROTEIN;
	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defo[num_proteins];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafo[num_proteins][num_proteins];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_pot[2][num_proteins];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_res[num_proteins];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_res[num_proteins];
	WorkerBuffer<dofMB<REAL>, HostPinnedAllocator<dofMB<REAL>>> h_dofSupport;
	WorkerBuffer<dofMB<REAL>, DeviceAllocator<dofMB<REAL>>> d_dofSupport;


	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	size_t blockSizeReduceRec = 0;
	size_t blockSizeReduceLig = 0;
	size_t dofSize[NUM_MAX_PROTEIN];
	ThreadSafeQueue<unsigned> stream_queue;
	StageResource _resource;
	std::mutex _lock;
	int numItemsInPipe; /** number of items in the pipeline */

	cudaStream_t stream; /** cuda streams */




};

template<typename REAL>
auto GPUEnergyServiceMBModes<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		//(void)numWorkers;
		std::vector<id_t> ids = { common->tableId};
		for ( int i = 0; i< common->numProteins; ++i){
			ids.push_back( common->proteins[i].gridId);
			ids.push_back( common->proteins[i].proteinId);
		}

		auto id = this->_dataMng->getCommonDeviceIds(ids);
		std::vector<as::workerId_t> vec(numWorkers);

		std::iota(vec.begin(), vec.end(), id[0]);
		return vec;
	};
	return fncObj;
}

template<typename REAL>
auto GPUEnergyServiceMBModes<REAL>::createItemProcessor() -> itemProcessor_t {

	std::shared_ptr<Private> p = std::make_shared<Private>();
	//deviceId_t deviceId = _workerId++;
	deviceId_t deviceId = 0;
	itemProcessor_t fncObj = [this, deviceId, p] (workItem_t* item) -> bool {

		/* Set the device to work with */
		cudaVerify(cudaSetDevice(deviceId));

				/* reset the predicates for the actual iteration*/
				//p->resetPrediacatesForIteration();
		p->configureDevice();
		if (item != nullptr) {

			auto dI = createStageResource(item, deviceId);

			const auto itemSize = item->size();
			assert(itemSize > 0);


			for (int id_protein = 0; id_protein < Common_MB_Modes::numProteins; ++id_protein){
				unsigned const numAtoms = dI.proteins[id_protein]->numAtoms;
				unsigned const& dofSize = 13 + Common_MB_Modes::numModes[id_protein];
				p->resizeBuffersIfRequired( itemSize, numAtoms, dofSize , id_protein);
			}
			p->addItemAndLigandSize(dI);

		} else {
			return false;
			//p->stagesMngt.rotate();
		}
		p->score();
		//	p->swapBuffers();

		return !(p->pipelineEmpty());

	};

	return fncObj;
}

} // namespace as

#endif

#endif /* SRC_SERVICE_GPUENERGYSERVICEMB_TPP_ */

//			 d_deform(
//					BLSZ_TRAFO,
//					gridSizeCenter,
//					stream,
//					*stageResc.proteins[id_centerProtein],
//					d_dof.get(0),
//					it->size(),
//					id_centerProtein,
//					stageResc.proteins[id_centerProtein]->xPosOrigin,
//					stageResc.proteins[id_centerProtein]->yPosOrigin,
//					stageResc.proteins[id_centerProtein]->zPosOrigin,
//					d_defo[0][id_centerProtein].getX(),
//					d_defo[0][id_centerProtein].getY(),
//					d_defo[0][id_centerProtein].getZ()
//					);

//					 d_transform(
//						BLSZ_TRAFO,
//						gridSizeCenter,
//						stream,
//						d_dof.get(0),
//						it->size(),
//						stageResc.proteins[id_centerProtein]->numAtoms,
//						*stageResc.proteins[id_centerProtein],
//						*stageResc.proteins[id_partnerProtein],
//						id_centerProtein,
//						id_partnerProtein,
////						stageResc.proteins[id_centerProtein]->pivot,
////						stageResc.proteins[id_partnerProtein]->pivot,
//						d_defo[id_centerProtein].getX(),
//						d_defo[id_centerProtein].getY(),
//						d_defo[id_centerProtein].getZ(),
//						d_trafo[id_centerProtein][id_partnerProtein].getX(),
//						d_trafo[id_centerProtein][id_partnerProtein].getY(),
//						d_trafo[id_centerProtein][id_partnerProtein].getZ()
//						);

//			 DEBUG
//			cudaDeviceSynchronize();
			//cudaDeviceSynchronize();
//									size_t bufferSize = d_dof[pipeIdxDof[1]].bufferSize();
//									WorkerBuffer<dof_t> h_dof(4,bufferSize);
//									size_t cpySize = h_dof.bufferSize()*sizeof(dof_t);
//
//									std::cout << "bufferSize: " << bufferSize << " cpySize: " << cpySize << std::endl;
//									cudaMemcpy(h_dof.get(0),d_dof[pipeIdxDof[0]].get(0), cpySize, cudaMemcpyDeviceToHost);
//									cudaMemcpy(h_dof.get(1),d_dof[pipeIdxDof[1]].get(0), cpySize, cudaMemcpyDeviceToHost);
//									cudaMemcpy(h_dof.get(2),d_dof[pipeIdxDof[2]].get(0), cpySize, cudaMemcpyDeviceToHost);
//									cudaMemcpy(h_dof.get(3),d_dof[pipeIdxDof[3]].get(0), cpySize, cudaMemcpyDeviceToHost);
//									std::cout << " stage one " <<std::endl;
//									for(size_t i = 0; i < bufferSize; ++i) {
//										std::cout << 0<< " " <<h_dof.get(0)[i]  << std::endl<< std::endl;
//									}
//									for(size_t i = 0; i < bufferSize; ++i) {
//										std::cout << 1<< " " <<h_dof.get(1)[i]  << std::endl<< std::endl;
//									}
//									for(size_t i = 0; i < bufferSize; ++i) {
//										std::cout << 2<< " " <<h_dof.get(2)[i]  << std::endl<< std::endl;
//									}
//									for(size_t i = 0; i < bufferSize; ++i) {
//										std::cout << 3<< " " <<h_dof.get(3)[i]  << std::endl<< std::endl ;
//									}
//									std::cout <<std::endl<< std::endl;
//				d_rotateForces(
//						BLSZ_INTRPL,
//						gridSizeCenter,
//						stream,
//						id_centerProtein,
//						d_pot[0][id_centerProtein].getX(),
//						d_pot[0][id_centerProtein].getY(),
//						d_pot[0][id_centerProtein].getZ(),
//						d_pot[0][id_centerProtein].getW(),
//						d_pot[1][id_centerProtein].getX(),
//						d_pot[1][id_centerProtein].getY(),
//						d_pot[1][id_centerProtein].getZ(),
//						d_pot[1][id_centerProtein].getW(),
//						d_dof.get(0),
//						stageResc.proteins[id_centerProtein]->numAtoms,
//						it->size()


//			d_DOFPos(
//				BLSZ_INTRPL,
//				gridSizeRec,
//				streams[id_stream],
//				*stageResc.rec,
//				d_dof[id_stream].get(0),
//				it->size(), 0,
//				d_defoRec[id_stream].getX(),
//				d_defoRec[id_stream].getY(),
//				d_defoRec[id_stream].getZ(),
//				d_trafoRec[id_stream].getX(),
//				d_trafoRec[id_stream].getY(),
//				d_trafoRec[id_stream].getZ()
//				);

//d_rotateForces(
					//						BLSZ_INTRPL,
					//						gridSizeCenter,
					//						stream,
					//						id_centerProtein,
					//						d_pot[0][id_centerProtein].getX(),
					//						d_pot[0][id_centerProtein].getY(),
					//						d_pot[0][id_centerProtein].getZ(),
					//						d_pot[0][id_centerProtein].getW(),
					//						d_pot[1][id_centerProtein].getX(),
					//						d_pot[1][id_centerProtein].getY(),
					//						d_pot[1][id_centerProtein].getZ(),
					//						d_pot[1][id_centerProtein].getW(),
					//						d_dof.get(0),
					//						stageResc.proteins[id_centerProtein]->numAtoms,
					//						it->size()eins[id_centerProtein]->numAtoms;
