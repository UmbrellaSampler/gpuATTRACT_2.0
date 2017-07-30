/*
 * CPUEnergy.h
 *
 *  Created on: Aug 9, 2016
 *      Author: uwe
 */

#ifndef SRC_CPU6DENERGYSERVICE_H_
#define SRC_CPU6DENERGYSERVICE_H_

#include "CPUService.h"
#include "publicTypes.h"
#include "Types_6D.h"
#include "nativeTypesWrapper.h"

namespace as {

class DataManager;

template<typename SERVICE>
class Configurator_6D;

template<typename REAL>
class CPU_6D_EnergyService : public CPUService<Types_6D<REAL>> {
public:
	using typename CPUService<Types_6D<REAL>>::service_t;
	using typename TypeWrapper<REAL>::real_t;
	//ToDo: refactor to input_t
	using dof_t = typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;

	/* serves as a builder class*/
//	using configurator_t = Configurator_6D<CPU_6D_EnergyService<REAL>>;
	using configurator_t = Configurator_6D<service_t>;

	/* need public for instantiate the server in configurator */

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;

	itemProcessor_t createItemProcessor() override;

	void setDataManager(std::shared_ptr<DataManager> mng) {
		_dataMng = mng;
	}
private:

	class Buffer;

	std::shared_ptr<DataManager> _dataMng;

};

}  // namespace as




#endif /* SRC_CPU6DENERGYSERVICE_H_ */
