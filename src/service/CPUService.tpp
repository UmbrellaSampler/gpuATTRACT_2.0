/*
 * CPUService.tpp
 *
 *  Created on: Jul 18, 2016
 *      Author: uwe
 */

#ifndef SRC_CPUSERVICE_TPP_
#define SRC_CPUSERVICE_TPP_

#include <vector>
#include <numeric>
#include "CPUService.h"
#include "Allocator.h"

namespace as {

template<typename GenericTypes>
auto CPUService<GenericTypes>::createDistributor() -> distributor_t {
	distributor_t fncObj = [] (common_t const* common, size_t numWorkers) {
		(void)common;
		std::vector<as::workerId_t> vec(numWorkers);
		std::iota(vec.begin(), vec.end(), 0);
		return vec;
	};

	return fncObj;

}

template<typename GenericTypes>
void CPUService<GenericTypes>::initAllocators() {
	this->setInputAllocator(std::make_shared<HostAllocator<input_t>>());
	this->setResultAllocator(std::make_shared<HostAllocator<result_t>>());
}

} // namespace as


#endif /* SRC_CPUSERVICE_TPP_ */
