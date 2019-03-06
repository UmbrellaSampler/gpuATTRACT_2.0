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
#include "scoring_kernel.h"
#include "macros.h"

#include "ThreadSafeQueue.h"

#include <iomanip>
#include <limits>
#include <mutex>

//#define DEBUG
#ifdef DEBUG
#include <fstream>
#include <string>
#include "debug_functions.h"
#endif

namespace as {

template<typename REAL>
GPUEnergyService6DModes<REAL>::GPUEnergyService6DModes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds) :
	GPUEnergyService<Types_6D_Modes<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds), _threadsPerDevice(1)
{}

template<typename REAL>
GPUEnergyService6DModes<REAL>::GPUEnergyService6DModes(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds, uint threadsPerDevice) :
	GPUEnergyService<Types_6D_Modes<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds), _threadsPerDevice(threadsPerDevice)
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

	Private(int deviceId) : numItemsInPipe(0) {
		for (unsigned i = 0; i<num_streams; ++i) {
			int i1d;
			cudaSetDevice(deviceId);
			cudaGetDevice(&i1d);
			CUDA_CHECK(cudaStreamCreate(&streams[i]));
			unsigned id_stream = i;
						stream_queue.push( id_stream );

		}
	}
	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBuffer( size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numAtomsLig, size_t const& dofSizeRec, size_t const& dofSizeLig, unsigned const id_stream ) {
		const size_t atomBufferSizeRec = numDOFs*numAtomsRec;
		const size_t atomBufferSizeLig = numDOFs*numAtomsLig;

		d_defoRec[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeRec));
		d_defoLig[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeLig));
		d_potRec[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSizeRec));
		d_potLig[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSizeLig));
		d_resRec[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSizeRec*numDOFs));
		d_resLig[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1, dofSizeLig*numDOFs));
		h_resRec[id_stream] = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSizeRec*numDOFs));
		h_resLig[id_stream] = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1, dofSizeLig*numDOFs));
		d_dof[id_stream]    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
		d_trafoRec[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeRec));
		d_trafoLig[id_stream] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSizeLig));
	}

	size_t bufferSize( unsigned const id_stream ) const {
		return d_trafoLig[id_stream].bufferSize(  );
	}

	void addItemAndLigandSize(StageResource const& resc, unsigned const id_stream) {
		_lock.lock();
		_resources[id_stream] = resc;
		++numItemsInPipe;
		_lock.unlock();
	}

	void resizeBuffersIfRequired(size_t const& numDOFs, size_t const& numAtomsRec, size_t const& numAtomsLig, size_t const& dofSizeRec, size_t const& dofSizeLig , unsigned const id_stream) {
		if (numDOFs*numAtomsLig > d_trafoLig[id_stream].bufferSize() || numDOFs*numAtomsRec > d_trafoRec[id_stream].bufferSize()) {
			allocateBuffer(numDOFs, numAtomsRec, numAtomsLig, dofSizeRec, dofSizeLig, id_stream );
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

	void enqueueStream( unsigned id_stream){
		stream_queue.push( id_stream );
	}

	unsigned  get_idStream(){
		return stream_queue.waitAndPop();
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

	void score(unsigned const id_stream ) {
			/* check if new item enters the pipeline */


		auto const& stageResc = _resources[id_stream];
		auto* const it = stageResc.item;

		cudaVerify(cudaMemcpyAsync(d_dof[id_stream].get(0), it->inputBuffer(),
				it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, streams[id_stream]));

		const unsigned numElRec = it->size()*stageResc.rec->numAtoms;
		const unsigned numElLig = it->size()*stageResc.lig->numAtoms;

		assert( numElRec <= d_trafoRec[id_stream].bufferSize() );
		assert( numElLig <= d_trafoLig[id_stream].bufferSize() );
		double radius_cutoff = it->common()->radius_cutoff;
		/* Perform cuda kernel calls */
		size_t gridSizeRec = ( numElRec + BLSZ_TRAFO - 1) / BLSZ_TRAFO;
		size_t gridSizeLig = ( numElLig + BLSZ_TRAFO - 1) / BLSZ_TRAFO;


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
//
//			d_DOFPos(
//				BLSZ_INTRPL,
//				gridSizeLig,
//				streams[id_stream],
//				*stageResc.lig,
//				d_dof[id_stream].get(0),
//				it->size(), 1,
//				d_defoLig[id_stream].getX(),
//				d_defoLig[id_stream].getY(),
//				d_defoLig[id_stream].getZ(),
//				d_trafoLig[id_stream].getX(),
//				d_trafoLig[id_stream].getY(),
//				d_trafoLig[id_stream].getZ()
//				);

//

			d_score(
				BLSZ_TRAFO,
				gridSizeRec,
				streams[id_stream],
				stageResc.gridLig.inner,
				stageResc.gridLig.outer,
				*stageResc.rec,
				d_dof[id_stream].get(0),
				it->size(),
				0,
				radius_cutoff,
				d_defoRec[id_stream].getX(),
				d_defoRec[id_stream].getY(),
				d_defoRec[id_stream].getZ(),
				d_trafoRec[id_stream].getX(),
				d_trafoRec[id_stream].getY(),
				d_trafoRec[id_stream].getZ(),
				d_potRec[id_stream].getX(),
				d_potRec[id_stream].getY(),
				d_potRec[id_stream].getZ(),
				d_potRec[id_stream].getW());

			d_score(
				BLSZ_TRAFO,
				gridSizeLig,
				streams[id_stream],
				stageResc.gridRec.inner,
				stageResc.gridRec.outer,
				*stageResc.lig,
				d_dof[id_stream].get(0),
				it->size(),
				1,
				radius_cutoff,
				d_defoLig[id_stream].getX(),
				d_defoLig[id_stream].getY(),
				d_defoLig[id_stream].getZ(),
				d_trafoLig[id_stream].getX(),
				d_trafoLig[id_stream].getY(),
				d_trafoLig[id_stream].getZ(),
				d_potLig[id_stream].getX(),
				d_potLig[id_stream].getY(),
				d_potLig[id_stream].getZ(),
				d_potLig[id_stream].getW());


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


#define DEBUG
#ifdef DEBUG1
			cudaDeviceSynchronize();
			std::string file_defo = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_defo += "/rec_defo.dat";
			std::fstream fs_defo;
			fs_defo.open (file_defo, std::fstream::in | std::fstream::out | std::fstream::trunc );
			size_t bufferSizeDefoRec1 = d_defoRec[id_stream].bufferSize();
			WorkerBuffer<REAL> h_DefoRec(3,bufferSizeDefoRec1);
			size_t cpySizeDefoRec1 = h_DefoRec.bufferSize()*sizeof(REAL);

			//std::cout << "bufferSize: " << bufferSizeDefoRec << " cpySize: " << cpySizeDefoRec << std::endl;
			cudaMemcpy(h_DefoRec.getX(),d_defoRec[id_stream].getX(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_DefoRec.getY(),d_defoRec[id_stream].getY(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_DefoRec.getZ(),d_defoRec[id_stream].getZ(), cpySizeDefoRec1, cudaMemcpyDeviceToHost);
			fs_defo << "x y z"<< std::endl;
			for(size_t i = 0; i < stageResc.rec->numAtoms; ++i) {
				fs_defo  << std::setprecision(10)<< h_DefoRec.getX()[i] << " " << h_DefoRec.getY()[i] << " " << h_DefoRec.getZ()[i] << std::endl;
			}
			fs_defo.close();

			cudaDeviceSynchronize();
			std::string file_trafor = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_trafor += "/rec_trafo.dat";
			std::fstream fs_trafor;
			fs_trafor.open (file_trafor, std::fstream::in | std::fstream::out | std::fstream::trunc );
			size_t bufferSizeTrafoRec1 = d_trafoRec[id_stream].bufferSize();
			WorkerBuffer<REAL> h_TrafoRec(3,bufferSizeTrafoRec1);
			size_t cpySizeTrafoRec1 = h_TrafoRec.bufferSize()*sizeof(REAL);

			//std::cout << "bufferSize: " << bufferSizeDefoRec << " cpySize: " << cpySizeDefoRec << std::endl;
			cudaMemcpy(h_TrafoRec.getX(),d_trafoRec[id_stream].getX(), cpySizeTrafoRec1, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_TrafoRec.getY(),d_trafoRec[id_stream].getY(), cpySizeTrafoRec1, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_TrafoRec.getZ(),d_trafoRec[id_stream].getZ(), cpySizeTrafoRec1, cudaMemcpyDeviceToHost);
			fs_trafor << "x y z"<< std::endl;
			for(size_t i = 0; i < stageResc.rec->numAtoms; ++i) {
				fs_trafor << std::setprecision(10) << std::setprecision(10)<< h_TrafoRec.getX()[i] << " " << h_TrafoRec.getY()[i] << " " << h_TrafoRec.getZ()[i] << std::endl;
			}
			fs_trafor.close();

			std::string file_trafo = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_trafo += "/lig_trafo.dat";
			std::fstream fs_trafo;
			fs_trafo.open (file_trafo, std::fstream::in | std::fstream::out | std::fstream::trunc );

			size_t bufferSizeTrafoLig = d_trafoLig[id_stream].bufferSize();

			WorkerBuffer<REAL> h_trafoLig(3,bufferSizeTrafoLig);
			size_t cpySizetrafoLig = h_trafoLig.bufferSize()*sizeof(REAL);
//
//			std::cout << "bufferSize: " << bufferSizeDefoLig << " cpySize: " << cpySizeDefoLig << std::endl;
			cudaMemcpy(h_trafoLig.getX(),d_trafoLig[id_stream].getX(), cpySizetrafoLig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_trafoLig.getY(),d_trafoLig[id_stream].getY(), cpySizetrafoLig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_trafoLig.getZ(),d_trafoLig[id_stream].getZ(), cpySizetrafoLig, cudaMemcpyDeviceToHost);
			fs_trafo << "x y z"<< std::endl;
			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
				fs_trafo << std::setprecision(10) << std::setprecision(10)<<h_trafoLig.getX()[i] << " " << h_trafoLig.getY()[i] << " " << h_trafoLig.getZ()[i] << std::endl;
			}
			fs_trafo.close();



			std::string fn_l_defo = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			fn_l_defo += "/lig_defo.dat";
			std::fstream f_l_defo;
			f_l_defo.open (fn_l_defo, std::fstream::in | std::fstream::out | std::fstream::trunc );

			size_t bufferSizeDefoLig = d_defoLig[id_stream].bufferSize();

			WorkerBuffer<REAL> h_defoLig(3,bufferSizeDefoLig);
			size_t cpySizedefoLig = h_defoLig.bufferSize()*sizeof(REAL);
//
//			std::cout << "bufferSize: " << bufferSizeDefoLig << " cpySize: " << cpySizeDefoLig << std::endl;
			cudaMemcpy(h_defoLig.getX(),d_defoLig[id_stream].getX(), cpySizedefoLig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_defoLig.getY(),d_defoLig[id_stream].getY(), cpySizedefoLig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_defoLig.getZ(),d_defoLig[id_stream].getZ(), cpySizedefoLig, cudaMemcpyDeviceToHost);
			f_l_defo << "x y z"<< std::endl;
			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
				f_l_defo << std::setprecision(10) <<h_defoLig.getX()[i] << " " << h_defoLig.getY()[i] << " " << h_defoLig.getZ()[i] << std::endl;
			}
			f_l_defo.close();




			std::string file_pot = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_pot += "/lig_pot.dat";
			std::fstream fs_pot;
			fs_pot.open (file_pot, std::fstream::in | std::fstream::out | std::fstream::trunc );
			WorkerBuffer<REAL> h_potLig(4,stageResc.lig->numAtoms);
			size_t cpySize = stageResc.lig->numAtoms*sizeof(REAL);
			//std::cout <<"fx fy fz"<<std::endl;
			fs_pot << "fx fy fz w"<< std::endl;
			cudaMemcpy(h_potLig.getX(),d_potLig[id_stream].getX(), cpySize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLig.getY(),d_potLig[id_stream].getY(), cpySize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLig.getZ(),d_potLig[id_stream].getZ(), cpySize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLig.getW(),d_potLig[id_stream].getW(), cpySize, cudaMemcpyDeviceToHost);
			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
				//			for(size_t i = 0; i < 20; ++i) {
				fs_pot<< std::setprecision(10) << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;//
			}
			fs_pot.close();


			std::string file_potr = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_potr += "/rec_pot.dat";
			std::fstream fs_potr;
			fs_potr.open (file_potr, std::fstream::in | std::fstream::out | std::fstream::trunc );
			WorkerBuffer<REAL> h_potRec(4,stageResc.rec->numAtoms);
			size_t cpySizerec = stageResc.rec->numAtoms*sizeof(REAL);
			//std::cout <<"fx fy fz"<<std::endl;
			fs_pot << "fx fy fz w"<< std::endl;
			cudaMemcpy(h_potRec.getX(),d_potRec[id_stream].getX(), cpySizerec, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potRec.getY(),d_potRec[id_stream].getY(), cpySizerec, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potRec.getZ(),d_potRec[id_stream].getZ(), cpySizerec, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potRec.getW(),d_potRec[id_stream].getW(), cpySizerec, cudaMemcpyDeviceToHost);
			for(size_t i = 0; i < stageResc.rec->numAtoms; ++i) {
				//			for(size_t i = 0; i < 20; ++i) {
				fs_potr<< std::setprecision(10) << h_potRec.getX()[i] << " " << h_potRec.getY()[i] << " " << h_potRec.getZ()[i] << " " << h_potRec.getW()[i] << std::endl;//
			}
			fs_potr.close();
			/* Perform cuda kernel calls */
			cudaDeviceSynchronize();
#endif
			gridSizeRec = ( numElRec + BLSZ_INTRPL - 1) / BLSZ_INTRPL;
			gridSizeLig = ( numElLig + BLSZ_INTRPL - 1) / BLSZ_INTRPL;


//			d_potForce (
//				BLSZ_INTRPL,
//				gridSizeLig,
//				streams[id_stream],
//				stageResc.gridRec.inner,
//				stageResc.gridRec.outer,
//				*stageResc.lig,
//				it->size(),
//				d_trafoLig[id_stream].getX(),
//				d_trafoLig[id_stream].getY(),
//				d_trafoLig[id_stream].getZ(),
//				d_potLig[id_stream].getX(),
//				d_potLig[id_stream].getY(),
//				d_potLig[id_stream].getZ(),
//				d_potLig[id_stream].getW()); // OK
//
//			d_potForce (
//				BLSZ_INTRPL,
//				gridSizeRec,
//				streams[id_stream],
//				stageResc.gridLig.inner,
//				stageResc.gridLig.outer,
//				*stageResc.rec,
//				it->size(),
//				d_trafoRec[id_stream].getX(),
//				d_trafoRec[id_stream].getY(),
//				d_trafoRec[id_stream].getZ(),
//				d_potRec[id_stream].getX(),
//				d_potRec[id_stream].getY(),
//				d_potRec[id_stream].getZ(),
//				d_potRec[id_stream].getW()); // OK

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
//			exit(EXIT_SUCCESS);
			d_NLPotForce<REAL, dof_t>(
				BLSZ_INTRPL,
				gridSizeLig,
				streams[id_stream],
				stageResc.gridRec.NL,
				*stageResc.rec,
				*stageResc.lig,
				*stageResc.table,
				radius_cutoff,
				*stageResc.simParam,
				it->size(),
				d_defoRec[id_stream].getX(),
				d_defoRec[id_stream].getY(),
				d_defoRec[id_stream].getZ(),
				d_trafoLig[id_stream].getX(),
				d_trafoLig[id_stream].getY(),
				d_trafoLig[id_stream].getZ(),
				d_potLig[id_stream].getX(),
				d_potLig[id_stream].getY(),
				d_potLig[id_stream].getZ(),
				d_potLig[id_stream].getW()); // OK


//
//			//get nl forces on receptor
			d_NLPotForce<REAL, dof_t>(
				BLSZ_INTRPL,
				gridSizeRec,
				streams[id_stream],
				stageResc.gridLig.NL,
				*stageResc.lig,
				*stageResc.rec,
				*stageResc.table,
				radius_cutoff,
				*stageResc.simParam,
				it->size(),
				d_defoLig[id_stream].getX(),
				d_defoLig[id_stream].getY(),
				d_defoLig[id_stream].getZ(),
				d_trafoRec[id_stream].getX(),
				d_trafoRec[id_stream].getY(),
				d_trafoRec[id_stream].getZ(),
				d_potRec[id_stream].getX(),
				d_potRec[id_stream].getY(),
				d_potRec[id_stream].getZ(),
				d_potRec[id_stream].getW()
				); // OK




			d_rotateForces(
				BLSZ_INTRPL,
				gridSizeRec,
				streams[id_stream],
				d_potRec[id_stream].getX(),
				d_potRec[id_stream].getY(),
				d_potRec[id_stream].getZ(),
				d_dof[id_stream].get(0),
				stageResc.rec->numAtoms,
				it->size()
				);
//			std::cout <<"after nl"<< std::endl;
#ifdef DEBUG1
			cudaDeviceSynchronize();
			std::string file_nl = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_nl += "/rec_nl.dat";
			std::fstream fs_nl;
			fs_nl.open (file_nl, std::fstream::in | std::fstream::out | std::fstream::trunc );
			cudaDeviceSynchronize();
			WorkerBuffer<REAL> h_potRecnl(4,stageResc.rec->numAtoms);
			size_t cpySizenlrec = stageResc.rec->numAtoms*sizeof(REAL);
			//std::cout <<"fx fy fz"<<std::endl;
			cudaMemcpy(h_potRecnl.getX(),d_potRec[id_stream].getX(), cpySizenlrec, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potRecnl.getY(),d_potRec[id_stream].getY(), cpySizenlrec, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potRecnl.getZ(),d_potRec[id_stream].getZ(), cpySizenlrec, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potRecnl.getW(),d_potRec[id_stream].getW(), cpySizenlrec, cudaMemcpyDeviceToHost);
			fs_nl << "fx fy fz w"<< std::endl;
			for(size_t i = 0; i < stageResc.rec->numAtoms; ++i) {
				//			for(size_t i = 0; i < 20; ++i) {
				fs_nl << std::setprecision(10)<< h_potRecnl.getX()[i] << " " << h_potRecnl.getY()[i] << " " << h_potRecnl.getZ()[i] << " " << h_potRecnl.getW()[i] << std::endl;

			}
			fs_nl.close();



			std::string file_nllig = getDebugPath<REAL, DOF_6D_Modes<REAL>>(true);
			file_nllig += "/lig_nl.dat";
			std::fstream fs_nllig;
			fs_nllig.open (file_nllig, std::fstream::in | std::fstream::out | std::fstream::trunc );
			cudaDeviceSynchronize();
			WorkerBuffer<REAL> h_potLignl(4,stageResc.rec->numAtoms);
			size_t cpySizenllig = stageResc.lig->numAtoms*sizeof(REAL);
			//std::cout <<"fx fy fz"<<std::endl;
			cudaMemcpy(h_potLignl.getX(),d_potLig[id_stream].getX(), cpySizenllig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLignl.getY(),d_potLig[id_stream].getY(), cpySizenllig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLignl.getZ(),d_potLig[id_stream].getZ(), cpySizenllig, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLignl.getW(),d_potLig[id_stream].getW(), cpySizenllig, cudaMemcpyDeviceToHost);
			fs_nllig << "fx fy fz w"<< std::endl;
			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
				//			for(size_t i = 0; i < 20; ++i) {
				fs_nllig << std::setprecision(10)<< h_potLignl.getX()[i] << " " << h_potLignl.getY()[i] << " " << h_potLignl.getZ()[i] << " " << h_potLignl.getW()[i] << std::endl;

			}
			fs_nllig.close();
#endif
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

			deviceReduce<REAL, DOF_6D_Modes<REAL>,0>(
				blockSizeReduceRec,
				it->size(),
				stageResc.rec,
				d_dof[id_stream].get(0),
				d_defoRec[id_stream].getX(), d_defoRec[id_stream].getY(), d_defoRec[id_stream].getZ(),
				d_potRec[id_stream].getX(), d_potRec[id_stream].getY(), d_potRec[id_stream].getZ(),
				d_potRec[id_stream].getW(),
				d_resRec[id_stream].get(0),
				streams[id_stream]);

			deviceReduce<REAL, DOF_6D_Modes<REAL>, 1>(
				blockSizeReduceLig,
				it->size(),
				stageResc.lig,
				d_dof[id_stream].get(0),
				d_defoLig[id_stream].getX(), 	d_defoLig[id_stream].getY(), d_defoLig[id_stream].getZ(),
				d_potLig[id_stream].getX(), 	d_potLig[id_stream].getY(), d_potLig[id_stream].getZ(),
				d_potLig[id_stream].getW(),
				d_resLig[id_stream].get(0),
				streams[id_stream]);


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
			cudaVerify(cudaMemcpyAsync( h_resRec[id_stream].get(0), d_resRec[id_stream].get(0), dofSizeRec*it->size()*sizeof(REAL),
					cudaMemcpyDeviceToHost, streams[id_stream]));
			cudaVerify(cudaMemcpyAsync( h_resLig[id_stream].get(0), d_resLig[id_stream].get(0), dofSizeLig*it->size()*sizeof(REAL),
					cudaMemcpyDeviceToHost, streams[id_stream]));
			cudaVerify(cudaStreamSynchronize(streams[id_stream]));


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

			nvtxRangePushA("Host");

			h_finalReduce< REAL, 1, true>(
				it->size(),
				stageResc.lig,
				it->inputBuffer(),
				h_resLig[id_stream].get(0),
				it->resultBuffer());
				nvtxRangePop();

			h_finalReduce< REAL, 0, true>(
				it->size(),
				stageResc.rec,
				it->inputBuffer(),
				h_resRec[id_stream].get(0),
				it->resultBuffer());
			nvtxRangePop();

			/* Signal that result is in buffer */
			it->setProcessed();

			/* signal that this stage was executed within the current iteration */


			/* signal that one item has been passed the last stage */
			signalItemPassedLastStage();
	}
	static unsigned constexpr num_streams = 1;
	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoRec[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoLig[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoRec[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potRec[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_resRec[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_resLig[num_streams];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_resRec[num_streams];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_resLig[num_streams];


	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	size_t blockSizeReduceRec = 0;
	size_t blockSizeReduceLig = 0;
	size_t dofSizeRec = 13 + Common_Modes::numModesRec;
	size_t dofSizeLig = 13 + Common_Modes::numModesLig;
	ThreadSafeQueue<unsigned> stream_queue;
	StageResource _resources[num_streams];
	std::mutex _lock;
	int numItemsInPipe; /** number of items in the pipeline */

	cudaStream_t streams[num_streams]; /** cuda streams */




};

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		//(void)numWorkers;
		std::vector<id_t> ids = {common->gridIdRec, common->gridIdLig, common->ligId, common->recId, common->tableId};
		auto id = this->_dataMng->getCommonDeviceIds(ids);
		std::vector<as::workerId_t> vec(numWorkers);
		std::iota(vec.begin(), vec.end(), 0);
		return vec;
	};
	return fncObj;
}

template<typename REAL>
auto GPUEnergyService6DModes<REAL>::createItemProcessor() -> itemProcessor_t {


	deviceId_t deviceId = _deviceIds[_workerId/_threadsPerDevice];
	std::shared_ptr<Private> p = std::make_shared<Private>(deviceId);
	 _workerId++;
	itemProcessor_t fncObj = [this, deviceId, p] (workItem_t* item) -> bool {

		/* Set the device to work with */
		cudaVerify(cudaSetDevice(deviceId));

				/* reset the predicates for the actual iteration*/
				//p->resetPrediacatesForIteration();
		p->configureDevice();
		unsigned id_stream = p->get_idStream();
		if (item != nullptr) {

			auto dI = createStageResource(item, deviceId);

			const auto itemSize = item->size();
			assert(itemSize > 0);

			const unsigned& numAtomsRec = dI.rec->numAtoms;
			unsigned const& numAtomsLig = dI.lig->numAtoms;
			unsigned const& dofSizeRec = 13 + Common_Modes::numModesRec;
			unsigned const& dofSizeLig = 13 + Common_Modes::numModesLig;

			p->addItemAndLigandSize(dI, id_stream);
			p->resizeBuffersIfRequired( itemSize, numAtomsRec, numAtomsLig, dofSizeRec, dofSizeLig , id_stream);



		} else {
			return false;
			//p->stagesMngt.rotate();
		}
		p->score( id_stream );
		p->enqueueStream( id_stream);
		//	p->swapBuffers();

		return !(p->pipelineEmpty());

	};

	return fncObj;
}

} // namespace as

#endif

#endif /* SRC_SERVICE_GPUENERGYSERVICE6D_TPP_ */
