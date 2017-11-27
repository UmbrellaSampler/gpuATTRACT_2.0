/*
 * GPUService.h
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_GPUSERVICE_H_
#define SRC_GPUSERVICE_H_

#ifdef CUDA

#include "EnergyService.h"

namespace as {

template<typename GenericTypes>
class GPUEnergyService : public EnergyService<GenericTypes> {
protected:
	using service_t = EnergyService<GenericTypes>;
public:

	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;

	explicit GPUEnergyService(std::shared_ptr<DataManager> dataMng);

	virtual ~GPUEnergyService() {};

	virtual itemProcessor_t createItemProcessor() = 0;

	distributor_t createDistributor() = 0;

};

}  // namespace as

#include "GPUEnergyService.tpp"

#endif


#endif /* SRC_GPUSERVICE_H_ */
