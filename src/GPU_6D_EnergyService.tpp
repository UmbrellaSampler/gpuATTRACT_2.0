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

#include "macros.h"

namespace as {

template<typename REAL>
GPU_6D_EnergyService<REAL>::GPU_6D_EnergyService() : _workerId(0) {
	_d = new Private();
}

template<typename REAL>
GPU_6D_EnergyService<REAL>::~GPU_6D_EnergyService() {
	delete _d;
}

template<typename REAL>
class GPU_6D_EnergyService<REAL>::Private {
	using dof_t = typename GPU_6D_EnergyService<REAL>::dof_t;
public:

	Private() : stagesMngt(numStages), LigMngt(numStages), pipeIdx{0,1}, numItemsInPipe(0) {
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
	void allocateBuffer(size_t numDOFs, size_t numAtoms) {
		const size_t atomBufferSize = numDOFs*numAtoms;
		for (int i = 0; i < 2; ++i) {
			d_dof[i]    = std::move(WorkerBuffer<dof_t,DeviceAllocator<dof_t>>(1,numDOFs));
			d_potLig[i] = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));
			d_res[i]    = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,13*atomBufferSize));
			h_res[i]    = std::move(WorkerBuffer<REAL, HostPinnedAllocator<REAL>>(4,13*atomBufferSize));
		}
		d_trafoLig = std::move(WorkerBuffer<REAL, DeviceAllocator<REAL>>(4,atomBufferSize));

	}

	WorkerBuffer<dof_t, DeviceAllocator<dof_t>> d_dof[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_trafoLig;
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_potLig[2];
	WorkerBuffer<REAL, DeviceAllocator<REAL>> d_res[2];
	WorkerBuffer<REAL, HostPinnedAllocator<REAL>> h_res[2];

	constexpr unsigned numStages = 5;
	bool predicates[2][numStages];
	RingArray<WorkerItem*> stagesMngt;
	RingArray<Protein*> LigMngt;
	unsigned pipeIdx[2];

	int numItemsInPipe; /** number of items in the pipeline */

	cudaStream_t streams[4]; /** cuda streams */

	cudaEvent_t events[7];   /** cuda events */

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

	deviceId_t deviceId = _workerId++;
	itemProcessor_t fncObj = [this, deviceId] (workItem_t* item) -> bool {

		/* Set the device to work with */
		CUDA_CHECK(cudaSetDevice(deviceId));

		assert(item->size() > 0);
		const auto itemSize = item->size();

		/* item pointers */
		const auto dofs = item->inputBuffer();
		const auto common = item->common();
		auto results = item->resultBuffer();

		/* get DataItem pointers */
		const auto grid = std::dynamic_pointer_cast<DeviceGridUnion<REAL>>(_dataMng->get(common->gridId, deviceId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(grid != nullptr);

		const auto lig = std::dynamic_pointer_cast<DeviceProtein<REAL>>(_dataMng->get(common->ligId, deviceId)).get();
		assert(lig != nullptr);

		const auto rec = std::dynamic_pointer_cast<DeviceProtein<REAL>>(_dataMng->get(common->recId, deviceId)).get();
		assert(rec != nullptr);

		const auto table = std::dynamic_pointer_cast<DeviceParamTable<REAL>>(_dataMng->get(common->tableId, deviceId)).get();
		assert(table != nullptr);

		const auto simParams = std::dynamic_pointer_cast<SimParam<REAL>>(_dataMng->get(common->paramsId, deviceId)).get();
		assert(simParams != nullptr);

	};

	return fncObj;
}

} // namespace as

#endif

#endif /* SRC_GPU_6D_ENERGYSERVICE_TPP_ */
