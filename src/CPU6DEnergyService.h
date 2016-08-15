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
#include "Types_2B_6D.h"
#include "nativeTypesWrapper.h"

namespace as {

class DataManager;

template<typename REAL>
class CPU_6D_EnergyService : public CPUService<typename Types_6D<REAL>::DOF, typename Types_6D<REAL>::Common,  typename Types_6D<REAL>::Result> {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using dof_t = typename Types_6D<real_t>::DOF;
	using common_t = typename Types_6D<real_t>::Common;
	using result_t = typename Types_6D<real_t>::Result;
	using service_t = CPUService<dof_t, common_t, result_t>;
	using typename service_t::workItem_t;

public:
	using typename service_t::itemProcessor_t;

	CPU_6D_EnergyService();
	CPU_6D_EnergyService(CPU_6D_EnergyService const&) = delete;
	CPU_6D_EnergyService& operator= (CPU_6D_EnergyService const&) = delete;
	CPU_6D_EnergyService(CPU_6D_EnergyService&& copy) {
		*this = std::move(copy);
	}
	CPU_6D_EnergyService& operator= (CPU_6D_EnergyService&& copy) {
		_d = std::move(copy._d);
		_dataMng = std::move(copy._dataMng);
		return *this;
	}
	virtual ~CPU_6D_EnergyService();

	itemProcessor_t createItemProcessor() override;

	void setDataManager(std::shared_ptr<DataManager> mng) {
		_dataMng = mng;
	}
private:

	class Private;
	Private* _d;
	std::shared_ptr<DataManager> _dataMng;

};

}  // namespace as




#endif /* SRC_CPU6DENERGYSERVICE_H_ */
