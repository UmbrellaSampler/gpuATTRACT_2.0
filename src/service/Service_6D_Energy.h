/*
 * Service_6D_Energy.h
 *
 *  Created on: Jul 30, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_SERVICE_6D_ENERGY_H_
#define SRC_SERVICE_SERVICE_6D_ENERGY_H_


#include "Service.h"
#include "Types_6D.h"

namespace as {

template<typename REAL>
class Service_6D_Energy : public Service<Types_6D<REAL>> {
protected:
	using service_t = Service<Types_6D<REAL>>;
public:

	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;

	virtual ~Service_6D_Energy() {};

	virtual itemProcessor_t createItemProcessor() = 0;

	distributor_t createDistributor() = 0;

	void initAllocators() = 0;

};


#endif /* SRC_SERVICE_SERVICE_6D_ENERGY_H_ */
