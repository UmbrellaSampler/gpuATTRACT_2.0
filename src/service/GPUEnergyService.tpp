/*
 * GPUService.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_GPUSERVICE_TPP_
#define SRC_GPUSERVICE_TPP_

#ifdef CUDA

#include <vector>
#include <numeric>
#include "GPUEnergyService.h"
#include "Allocator.h"

namespace as {

template<typename GenericTypes>
GPUEnergyService<GenericTypes>::GPUEnergyService(std::shared_ptr<DataManager> dataMng) :
	EnergyService<GenericTypes>(dataMng)
{}

template<typename GenericTypes>
void GPUEnergyService<GenericTypes>::initAllocators() {
	this->setInputAllocator(std::make_shared<HostPinnedAllocator<input_t>>());
	this->setResultAllocator(std::make_shared<HostPinnedAllocator<result_t>>());
}

} // namespace as

#endif



#endif /* SRC_GPUSERVICE_TPP_ */
