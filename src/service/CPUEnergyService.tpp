/*
 * CPUService.tpp
 *
 *  Created on: Jul 18, 2016
 *      Author: uwe
 */

#ifndef SRC_CPUSERVICE_TPP_
#define SRC_CPUSERVICE_TPP_

#include "CPUEnergyService.h"

#include <vector>
#include <numeric>
#include "Allocator.h"

namespace as {

template<typename GenericTypes>
CPUEnergyService<GenericTypes>::CPUEnergyService(std::shared_ptr<DataManager> dataMng) :
	EnergyService<GenericTypes>(std::make_shared<HostAllocator<input_t>>(),
			std::make_shared<HostAllocator<result_t>>(), dataMng)
{}

template<typename GenericTypes>
auto CPUEnergyService<GenericTypes>::createDistributor() -> distributor_t {
	distributor_t fncObj = [] (common_t const* common, std::size_t numWorkers) {
		(void)common;
		std::vector<as::workerId_t> vec(numWorkers);
		std::iota(vec.begin(), vec.end(), 0);
		return vec;
	};

	return fncObj;

}

} // namespace as


#endif /* SRC_CPUSERVICE_TPP_ */
