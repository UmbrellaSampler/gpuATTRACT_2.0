/*
 * CPUEnergy.h
 *
 *  Created on: Aug 9, 2016
 *      Author: uwe
 */

#ifndef SRC_CPUENERGYSERVICE6D_H_
#define SRC_CPUENERGYSERVICE6D_H_

#include "CPUEnergyService.h"
#include "publicTypes.h"
#include "Types_6D.h"
#include "nativeTypesWrapper.h"

namespace as {

template<typename SERVICE>
class Configurator_6D;

template<typename REAL>
class CPUEnergyService6D : public CPUEnergyService<Types_6D<REAL>> {
public:
	using typename CPUEnergyService<Types_6D<REAL>>::service_t;
	using real_t = typename TypeWrapper<REAL>::real_t;
	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;
	using configurator_t = Configurator_6D<real_t>;

	/* need public for instantiate the server in configurator */

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;

	explicit CPUEnergyService6D(std::shared_ptr<DataManager> dataMng);
	virtual ~CPUEnergyService6D() {};

	itemProcessor_t createItemProcessor() override;

private:

	class Buffer;

};

}  // namespace as




#endif /* SRC_CPUENERGYSERVICE6D_H_ */
