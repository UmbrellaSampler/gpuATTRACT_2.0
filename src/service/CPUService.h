/*
 * CPUService.h
 *
 *  Created on: Jul 18, 2016
 *      Author: uwe
 */

#ifndef SRC_CPUSERVICE_H_
#define SRC_CPUSERVICE_H_

#include "Service.h"

namespace as {

//template<typename InputType, typename CommonType, typename ResultType>
//class CPUService : public Service<InputType,CommonType,ResultType> {
template<typename GenericTypes>
class CPUService : public Service<GenericTypes> {
protected:
	using service_t = Service<GenericTypes>;
public:

	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;

	virtual ~CPUService() {};

	virtual itemProcessor_t createItemProcessor() = 0;

	distributor_t createDistributor() override;

	void initAllocators() override;

};

}  // namespace as

#include "CPUService.tpp"


#endif /* SRC_CPUSERVICE_H_ */
