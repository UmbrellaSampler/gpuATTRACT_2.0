/*
 * GPU_6D_EnergyService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_GPU_6D_ENERGYSERVICE_TPP_
#define SRC_GPU_6D_ENERGYSERVICE_TPP_

#ifdef CUDA

#include "GPU_6D_EnergyService.h"

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

#include "macros.h"

// ToDo: remove
#include <iostream>

namespace as {

template<typename REAL>
struct StageResource {
private:
	using workItem_t = typename GPU_6D_EnergyService<REAL>::workItem_t;
public:
	d_GridUnion<REAL> grid;
	d_Protein<REAL>* rec;
	d_Protein<REAL>* lig;
	d_ParamTable<REAL>* table;
	SimParam<REAL>* simParam;
	workItem_t* item;
};

template<typename REAL>
class GPU_6D_EnergyService<REAL>::Private {
	using dof_t = typename GPU_6D_EnergyService<REAL>::dof_t;
	using workItem_t = typename GPU_6D_EnergyService<REAL>::workItem_t;
public:

	Private() : stagesMngt(numStages), pipeIdx{0,1}, numItemsInPipe(0) {
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
	void allocateBuffer(size_t const& numDOFs, size_t const& numAtoms) {
		const size_t atomBufferSize = numDOFs*numAtoms;
		for (int i = 0; i < 2; ++i) {
			d_dof[i]    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
			d_potLig[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
			d_res[i]    = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(1,13*numDOFs));
			h_res[i]    = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(1,13*numDOFs));
		}
		d_trafoLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(3,atomBufferSize));
	}

	size_t bufferSize() {
		return d_trafoLig.bufferSize();
	}

	void addItemAndLigandSize(StageResource<REAL> const& resc) {
		++numItemsInPipe;
		predicates[pipeIdx[0]][0] = true;
		stagesMngt.push(resc);
	}

	void resizeBuffersIfRequired(size_t const& numDOFs, size_t const& numAtoms) {
		if (numDOFs*numAtoms > bufferSize()) {
			allocateBuffer(numDOFs, numAtoms);
		}
	}

	void swapBuffers() {
		std::swap(pipeIdx[0], pipeIdx[1]);
	}

	bool pipelineEmpty() {
		return numItemsInPipe == 0;
	}

	void S1_copyH2D() {
		if (predicates[pipeIdx[0]][0])
		{
			constexpr unsigned stageId = 0;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

			cudaVerify(cudaStreamWaitEvent(streams[0], events[2], 0));
			cudaVerify(cudaMemcpyAsync(d_dof[pipeIdx[0]].get(0), it->inputBuffer(),
					it->size()*sizeof(dof_t), cudaMemcpyHostToDevice, streams[0]));
			cudaVerify(cudaEventRecord(events[0], streams[0]));

		}
	}

	void S2_transform_potForce() {
		/* check if stage 0 was executed in last iteration */
		if (predicates[pipeIdx[1]][0])
		{
			constexpr unsigned stageId = 1;
			auto const& stageResc = stagesMngt.get(stageId);
			auto* const it = stageResc.item;

			const unsigned numEl = it->size()*stageResc.lig->numAtoms;
//			std::cout << "GPUWorker.cpp: " << "numEl " << numEl << " _atomBufferSize " << _atomBufferSize << std::endl;
			assert(numEl <= bufferSize());

			/* Perform cuda kernel calls */
			size_t gridSize = ( numEl + BLSZ_TRAFO - 1) / BLSZ_TRAFO;


			d_DOF2Pos(
					BLSZ_TRAFO,
					gridSize,
					streams[2],
					stageResc.lig->xPos,
					stageResc.lig->yPos,
					stageResc.lig->zPos,
					d_dof[pipeIdx[1]].getX(),
					stageResc.lig->numAtoms,
					it->size(),
					d_trafoLig.getX(),
					d_trafoLig.getY(),
					d_trafoLig.getZ());

			//DEBUG
//			cudaDeviceSynchronize();
//			WorkerBuffer<REAL> h_trafoLig(3,stageResc.lig->numAtoms);
//			size_t cpySize = stageResc.lig->numAtoms*sizeof(REAL);
//			cudaMemcpy(h_trafoLig.getX(),d_trafoLig.getX(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_trafoLig.getY(),d_trafoLig.getY(), cpySize, cudaMemcpyDeviceToHost);
//			cudaMemcpy(h_trafoLig.getZ(),d_trafoLig.getZ(), cpySize, cudaMemcpyDeviceToHost);
//			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
//				std::cout << h_trafoLig.getX()[i] << " " << h_trafoLig.getY()[i] << " " << h_trafoLig.getZ()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);



			/* Device: Wait for completion of copyH2D of DOFs to complete */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[0], 0));
			/* Device: Signal event when transformation has completed */
			cudaVerify(cudaEventRecord(events[2], streams[2]));
			/* Device: Wait for completion of reduction of the previous round */
			cudaVerify(cudaStreamWaitEvent(streams[2], events[5+pipeIdx[1]], 0));

			/* Perform cuda kernel calls */
			gridSize = ( numEl + BLSZ_INTRPL - 1) / BLSZ_INTRPL;
			d_potForce (
				BLSZ_INTRPL,
				gridSize,
				streams[2],
				stageResc.grid.inner,
				stageResc.grid.outer,
				*stageResc.lig,
				it->size(),
				d_trafoLig.getX(),
				d_trafoLig.getY(),
				d_trafoLig.getZ(),
				d_potLig[pipeIdx[1]].getX(),
				d_potLig[pipeIdx[1]].getY(),
				d_potLig[pipeIdx[1]].getZ(),
				d_potLig[pipeIdx[1]].getW());


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
				streams[2],
				stageResc.grid.NL,
				*stageResc.rec,
				*stageResc.lig,
				*stageResc.table,
				*stageResc.simParam,
				it->size(),
				d_trafoLig.getX(),
				d_trafoLig.getY(),
				d_trafoLig.getZ(),
				d_potLig[pipeIdx[1]].getX(),
				d_potLig[pipeIdx[1]].getY(),
				d_potLig[pipeIdx[1]].getZ(),
				d_potLig[pipeIdx[1]].getW());

			cudaDeviceSynchronize();
//			exit(EXIT_SUCCESS);
			WorkerBuffer<REAL> h_potLig(4,stageResc.lig->numAtoms);
			size_t cpySize = stageResc.lig->numAtoms*sizeof(REAL);
			cudaMemcpy(h_potLig.getX(),d_potLig[pipeIdx[1]].getX(), cpySize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLig.getY(),d_potLig[pipeIdx[1]].getY(), cpySize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLig.getZ(),d_potLig[pipeIdx[1]].getZ(), cpySize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_potLig.getW(),d_potLig[pipeIdx[1]].getW(), cpySize, cudaMemcpyDeviceToHost);
			for(size_t i = 0; i < stageResc.lig->numAtoms; ++i) {
//			for(size_t i = 0; i < 20; ++i) {
				std::cout << h_potLig.getX()[i] << " " << h_potLig.getY()[i] << " " << h_potLig.getZ()[i] << " " << h_potLig.getW()[i] << std::endl;
			}
			exit(EXIT_SUCCESS);

			/* Device: Signal event when force and energy calc. has completed */
			cudaVerify(cudaEventRecord(events[3+pipeIdx[1]], streams[2]));

			/* signal that this stage was executed within the current iteration */
			predicates[pipeIdx[0]][1] = true;
		}
	}

	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_res[2];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_res[2];

	static constexpr unsigned numStages = 5;
	bool predicates[2][numStages];
	RingArray<StageResource<REAL>> stagesMngt;
	unsigned pipeIdx[2];

	int numItemsInPipe; /** number of items in the pipeline */

	cudaStream_t streams[4]; /** cuda streams */

	cudaEvent_t events[7];   /** cuda events */

	static constexpr size_t BLSZ_TRAFO = 128;
	static constexpr size_t BLSZ_INTRPL = 128;
	static constexpr size_t BLSZ_REDUCE = 512;

};

template<typename REAL>
auto GPU_6D_EnergyService<REAL>::createDistributor() -> distributor_t {
	distributor_t fncObj = [this] (common_t const* common, size_t numWorkers) {
		(void)numWorkers;
		std::vector<id_t> ids = {common->gridId, common->ligId, common->recId, common->tableId};
		auto vec = _dataMng->getCommonDeviceIds(ids);
		return vec;
	};
	return fncObj;
}

template<typename REAL>
auto GPU_6D_EnergyService<REAL>::createItemProcessor() -> itemProcessor_t {

	std::shared_ptr<Private> p = std::make_shared<Private>();
	deviceId_t deviceId = _workerId++;
	itemProcessor_t fncObj = [this, deviceId, p] (workItem_t* item) -> bool {

		/* Set the device to work with */
		cudaVerify(cudaSetDevice(deviceId));

		/* reset the predicates for the actual iteration*/
		for (unsigned i = 0; i<p->numStages; ++i) {
			p->predicates[p->pipeIdx[0]][i] = false;
		}

		if (item != nullptr) {

			/* item pointers */
//			const auto dofs = item->inputBuffer();
			const auto common = item->common();
//			auto results = item->resultBuffer();

			/* get DataItem pointers */
			auto grid = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(_dataMng->get(common->gridId, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
			assert(grid != nullptr);

			auto lig = std::dynamic_pointer_cast<DeviceProtein<REAL>>(_dataMng->get(common->ligId, deviceId)).get();
			assert(lig != nullptr);

			auto rec = std::dynamic_pointer_cast<DeviceProtein<REAL>>(_dataMng->get(common->recId, deviceId)).get();
			assert(rec != nullptr);

			auto table = std::dynamic_pointer_cast<DeviceParamTable<REAL>>(_dataMng->get(common->tableId, deviceId)).get();
			assert(table != nullptr);

			auto simParam = std::dynamic_pointer_cast<SimParam<REAL>>(_dataMng->get(common->paramsId)).get();
			assert(simParam != nullptr);


			StageResource<REAL> dI;

			dI.grid 	= grid->getDesc();
			dI.lig 		= &lig->desc;
			dI.rec 		= &rec->desc;
			dI.table 	= &table->desc;
			dI.simParam = simParam;
			dI.item 	= item;

			const auto itemSize = item->size();
			assert(itemSize > 0);

			const unsigned& numAtoms = dI.lig->numAtoms;
			p->addItemAndLigandSize(dI);
			p->resizeBuffersIfRequired(itemSize, numAtoms);

		} else {
			p->stagesMngt.rotate();
		}

		p->S1_copyH2D();
		p->S2_transform_potForce();
		p->swapBuffers();

		return !(p->pipelineEmpty());

	};

	return fncObj;
}

} // namespace as

#endif

#endif /* SRC_GPU_6D_ENERGYSERVICE_TPP_ */
