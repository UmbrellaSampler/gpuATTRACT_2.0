/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICE6D_TPP_
#define SRC_SERVICE_GPUENERGYSERVICE6D_TPP_

#ifdef CUDA

#include <nvToolsExt.h>
#include "GPUEnergyService6D.h"

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

#include "transform.h"
#include "interpolation.h"
#include "neighborlist.h"
#include "reduction.h"

#include "macros.h"

// ToDo: remove
#include <iostream>
#include <mutex>
#include "ThreadSafeQueue.h"
#include "scoring_kernel.h"
//#include "cudaProfiler.h"

namespace as {

template<typename REAL>
GPUEnergyService6D<REAL>::GPUEnergyService6D(std::shared_ptr<DataManager> dataMng,
		std::vector<int> const& deviceIds) :
	GPUEnergyService<Types_6D<REAL>>::GPUEnergyService(dataMng), _workerId(0), _deviceIds(deviceIds)
{}

template<typename REAL>
struct GPUEnergyService6D<REAL>::StageResource {
private:
	using workItem_t = typename GPUEnergyService6D<REAL>::workItem_t;
public:
	d_GridUnion<REAL> grid;
	d_Protein<REAL>* rec;
	d_Protein<REAL>* lig;
	d_ParamTable<REAL>* table;
	SimParam<REAL>* simParam;
	workItem_t* item;
};

template<typename REAL>
auto GPUEnergyService6D<REAL>::createStageResource(workItem_t* item, unsigned const& deviceId) -> StageResource {
	/* item pointers */
//			const auto dofs = item->inputBuffer();
	const auto common = item->common();
//			auto results = item->resultBuffer();

	/* get DataItem pointers */
	auto grid = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(this->_dataMng->get(common->gridId, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
	assert(grid != nullptr);

	auto lig = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->ligId, deviceId)).get();
	assert(lig != nullptr);

	auto rec = std::dynamic_pointer_cast<DeviceProtein<REAL>>(this->_dataMng->get(common->recId, deviceId)).get();
	assert(rec != nullptr);

	auto table = std::dynamic_pointer_cast<DeviceParamTable<REAL>>(this->_dataMng->get(common->tableId, deviceId)).get();
	assert(table != nullptr);

	auto simParam = std::dynamic_pointer_cast<SimParam<REAL>>(this->_dataMng->get(common->paramsId)).get();
	assert(simParam != nullptr);

	StageResource stageResource;
	stageResource.grid 	= grid->getDesc();
	stageResource.lig 		= &lig->desc;
	stageResource.rec 		= &rec->desc;
	stageResource.table 	= &table->desc;
	stageResource.simParam = simParam;
	stageResource.item 	= item;

	return stageResource;
}

template<typename REAL>
class GPUEnergyService6D<REAL>::Private {
	using dof_t = typename GPUEnergyService6D<REAL>::input_t;
	using workItem_t = typename GPUEnergyService6D<REAL>::workItem_t;

public:

	Private() :   numItemsInPipe(0) {
		for(unsigned i = 0; i< num_streams; ++i){
			unsigned id_stream = i;
			stream_queue.push( id_stream );
			CUDA_CHECK(cudaStreamCreate( &streams[i]) );
			CUDA_CHECK(cudaEventCreate( &events[i], cudaEventDisableTiming) );
		}
	}
	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBuffer(size_t const& numDOFs, size_t const& numAtoms, unsigned const id_buffer) {
		const size_t atomBufferSize = numDOFs*numAtoms;
		d_dof[id_buffer]    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
		d_potLig[id_buffer] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
		d_res[id_buffer]    = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1,13*numDOFs));
		h_res[id_buffer]    = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1,13*numDOFs));
		d_trafoLig[id_buffer]  = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSize));
	}

	size_t bufferSize(unsigned const id_stream ) const {
		return d_trafoLig[id_stream].bufferSize();
	}

	void addItemAndLigandSize(StageResource const& resc, unsigned const id_stream) {
		//_lock.lock();
		++numItemsInPipe;
		//_lock.unlock();
		//predicates[pipeIdx[0]][0] = true;
		resources[id_stream] = resc;
		//stagesMngt.push(resc);
	}

	void resizeBuffersIfRequired(size_t const& numDOFs, size_t const& numAtoms, unsigned const id_stream ) {
		if (numDOFs*numAtoms > bufferSize( id_stream )) {
			allocateBuffer(numDOFs, numAtoms, id_stream);
		}
	}

	bool pipelineEmpty()  {
		//_lock.lock();
		bool tmp = numItemsInPipe;
		//_lock.unlock();
		return tmp == 0;

	}

	void signalItemPassedLastStage() {
		//_lock.lock();
		--numItemsInPipe;
		//_lock.unlock();
	}

	void enqueueStream( unsigned id_stream){
		stream_queue.push( id_stream );
	}

	unsigned  get_idStream(){
		return stream_queue.waitAndPop();
	}

	void configureDevice() {

		if (blockSizeReduce == 0) {
			int id;
			cudaVerify(cudaGetDevice(&id));
			cudaDeviceProp deviceProp;
			cudaVerify(cudaGetDeviceProperties(&deviceProp, id));
			size_t sharedMem = deviceProp.sharedMemPerBlock;
			size_t pow2 = 16;
			while (pow2*13*sizeof(REAL) < sharedMem) {
				pow2 *= 2;
			}

			blockSizeReduce = pow2 / 2;
		}

	}
	int getnumpipe(){
		return numItemsInPipe;
	}

	void score(unsigned const id_stream) {
		/* check if new item enters the pipeline */

			auto const& stageResc = resources[id_stream];
			auto* const it = stageResc.item;
			const unsigned numEl = it->size()*stageResc.lig->numAtoms;
			assert(numEl <= bufferSize( id_stream ));


//			//DEBUG
//			for (size_t i = 0; i < it->size(); ++i) {
//				std::cout << it->inputBuffer()[i] << std::endl;
//			}

			cudaVerify(cudaMemcpyAsync(d_dof[id_stream].get(0), it->inputBuffer(),
					it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, streams[id_stream]));

			/* Perform cuda kernel calls */
			size_t gridSize = ( numEl + BLSZ_TRAFO - 1) / BLSZ_TRAFO;

			/* Device: Wait for completion of copyH2D of DOFs to complete */
			d_score(
				BLSZ_TRAFO,
				gridSize,
				streams[id_stream],
				stageResc.grid.inner,
				stageResc.grid.outer,
				*stageResc.lig,
				d_dof[id_stream].get(0),
				it->size(),
				1,
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

//			d_DOF2Pos(
//					BLSZ_TRAFO,
//					gridSize,
//					streams[id_stream],
//					stageResc.lig->xPos,
//					stageResc.lig->yPos,
//					stageResc.lig->zPos,
//					d_dof[id_stream].get(0),
//					stageResc.lig->numAtoms,
//					it->size(),
//					d_trafoLig[id_stream].getX(),
//					d_trafoLig[id_stream].getY(),
//					d_trafoLig[id_stream].getZ()); //OK

			// DEBUG
//			cudaDeviceSynchronize();
//			size_t bufferSize = d_trafoLig.bufferSize();
//			WorkerBuffer<REAL> h_trafoLig(3,bufferSize);
//			size_t cpySize = h_trafoLig.bufferSize()*sizeof(REAL);
//			std::cout << "#gpu trafo" << std::endl;
//			//std::cout << "bufferSize: " << bufferSize << " cpySize: " << cpySize << std::endl;
//			cudaMemcpy(h_trafoLig.getX(),d_trafoLig.getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_trafoLig.getY(),d_trafoLig.getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_trafoLig.getZ(),d_trafoLig.getZ(), cpySize, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < bufferSize; ++i) {
//				std::cout << h_trafoLig.getX()[i] << " " << h_trafoLig.getY()[i] << " " << h_trafoLig.getZ()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);
			/* Perform cuda kernel calls */

			gridSize = ( numEl + BLSZ_INTRPL - 1) / BLSZ_INTRPL;

//			d_potForce (
//				BLSZ_INTRPL,
//				gridSize,
//				streams[id_stream],
//				stageResc.grid.inner,
//				stageResc.grid.outer,
//				*stageResc.lig,
//				it->size(),
//				d_trafoLig[id_stream].getX(),
//				d_trafoLig[id_stream].getY(),
//				d_trafoLig[id_stream].getZ(),
//				d_potLig[id_stream].getX(),
//				d_potLig[id_stream].getY(),
//				d_potLig[id_stream].getZ(),
//				d_potLig[id_stream].getW()); // OK


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
				gridSize,
				streams[id_stream],
				stageResc.grid.NL,
				*stageResc.rec,
				*stageResc.lig,
				*stageResc.table,
				*stageResc.simParam,
				it->size(),
				d_trafoLig[id_stream].getX(),
				d_trafoLig[id_stream].getY(),
				d_trafoLig[id_stream].getZ(),
				d_potLig[id_stream].getX(),
				d_potLig[id_stream].getY(),
				d_potLig[id_stream].getZ(),
				d_potLig[id_stream].getW()); // OK

//			// Debug
//			cudaDeviceSynchronize();
//			size_t bufferSize = d_potLig[pipeIdx[1]].bufferSize();
//			WorkerBuffer<REAL> h_potLig(4,bufferSize);
//			size_t cpySize = h_potLig.bufferSize()*sizeof(REAL);
//			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[1]].getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[1]].getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[1]].getW(), cpySize, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < bufferSize; ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);


			/* Device: Wait for completion of PotForce calc. to complete */


			deviceReduce(
				blockSizeReduce,
				it->size(),
				stageResc.lig->numAtoms,
				stageResc.lig->xPos,
				stageResc.lig->yPos,
				stageResc.lig->zPos,
				d_potLig[id_stream].getX(),
				d_potLig[id_stream].getY(),
				d_potLig[id_stream].getZ(),
				d_potLig[id_stream].getW(),
				d_res[id_stream].get(0),
				streams[id_stream]);

//			cudaDeviceSynchronize();
//			unsigned numDofs = it->size();
//			WorkerBuffer<REAL> h_potLig(1,13*numDofs);
//			size_t cpySize = h_potLig.bufferSize()*sizeof(REAL);
//			cudaMemcpy(h_potLig.get(0),d_res[pipeIdx[0]].get(0), cpySize, cudaMemcpyDeviceToHost);
//
//			for(size_t i = 0; i < numDofs; ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				REAL x = h_potLig.get(0)[i*13 + 0];
//				REAL y = h_potLig.get(0)[i*13 + 1];
//				REAL z = h_potLig.get(0)[i*13 + 2];
//				REAL E = h_potLig.get(0)[i*13 + 3];
//				std::cout << x << " " << y << " " << z << " " << E << std::endl;
//			}
//			exit(EXIT_SUCCESS);


			/* copy results to host */
			cudaVerify(cudaMemcpyAsync(h_res[id_stream].get(0), d_res[id_stream].get(0), 13*it->size()*sizeof(REAL),
					cudaMemcpyDeviceToHost, streams[id_stream]));


			/* Device: Signal event when reduction has completed */

//			cudaVerify(cudaEventRecord(events[id_stream], streams[id_stream]));
//			cudaVerify(cudaEventSynchronize(events[id_stream]));
			cudaVerify(cudaStreamSynchronize(streams[id_stream]));

			nvtxRangePushA("Host");
			h_finalReduce(
				it->size(),
				it->inputBuffer(),
				h_res[id_stream].get(0),
				it->resultBuffer());
			nvtxRangePop();

			/* Signal that result is in buffer */
			it->setProcessed();
			/* signal that one item has been passed the last stage */
			signalItemPassedLastStage();
			/* signal that this stage was executed within the current iteration */


	}


	static unsigned constexpr num_streams = 1;
	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_defoLig[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig[num_streams];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_res[num_streams];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_res[num_streams];

	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	size_t blockSizeReduce = 0;
	static constexpr unsigned numStages = 5;


	int numItemsInPipe; /** number of items in the pipeline */
	StageResource resources[num_streams];
	cudaStream_t streams[num_streams]; /** cuda streams */
	ThreadSafeQueue<unsigned> stream_queue;
	cudaEvent_t events[7];   /** cuda events */
	std::mutex _lock;


};

template<typename REAL>
auto GPUEnergyService6D<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		//(void)numWorkers;
		std::vector<id_t> ids = {common->gridId, common->ligId, common->recId, common->tableId};

		auto id = this->_dataMng->getCommonDeviceIds(ids);
		std::vector<as::workerId_t> vec(numWorkers);

		std::iota(vec.begin(), vec.end(), id[0]);
		return vec;
	};
	return fncObj;
}

template<typename REAL>
auto GPUEnergyService6D<REAL>::createItemProcessor() -> itemProcessor_t {

	std::shared_ptr<Private> p = std::make_shared<Private>();
	//deviceId_t deviceId = _workerId++;
	deviceId_t deviceId = 0;
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

			const unsigned& numAtoms = dI.lig->numAtoms;
			p->addItemAndLigandSize(dI, id_stream);
			p->resizeBuffersIfRequired(itemSize, numAtoms,id_stream);

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
