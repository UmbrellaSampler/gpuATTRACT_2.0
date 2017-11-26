/*
 * CPUService.h
 *
 *  Created on: Jul 18, 2016
 *      Author: uwe
 */

#ifndef SRC_CPUSERVICE_H_
#define SRC_CPUSERVICE_H_

#include "EnergyService.h"

namespace as {

//template<typename InputType, typename CommonType, typename ResultType>
//class CPUService : public Service<InputType,CommonType,ResultType> {
template<typename GenericTypes>
class CPUEnergyService : public EnergyService<GenericTypes> {
protected:
public:
	using typename EnergyService<GenericTypes>::service_t;

	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;

	explicit CPUEnergyService(std::shared_ptr<DataManager> dataMng);
	virtual ~CPUEnergyService() {};

	virtual itemProcessor_t createItemProcessor() = 0;

	distributor_t createDistributor() override;

};

}  // namespace as

#include "CPUEnergyService.tpp"


#endif /* SRC_CPUSERVICE_H_ */
