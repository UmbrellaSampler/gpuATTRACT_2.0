/*
 * CPUEnergy.h
 *
 *  Created on: Aug 9, 2016
 *      Author: uwe
 */

#ifndef SRC_CPUENERGYSERVICE_H_
#define SRC_CPUENERGYSERVICE_H_

#include "CPUService.h"
#include "publicTypes.h"
#include "Types_2B_6D.h"
#include "nativeTypesWrapper.h"

namespace as {

class DataManager;

template<typename REAL>
class CPUEnergyService : public CPUService<typename Types_2B_6D<REAL>::DOF, typename Types_2B_6D<REAL>::Common,  typename Types_2B_6D<REAL>::Result> {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using dof_t = typename Types_2B_6D<real_t>::DOF;
	using common_t = typename Types_2B_6D<real_t>::Common;
	using result_t = typename Types_2B_6D<real_t>::Result;
	using service_t = CPUService<dof_t, common_t, result_t>;
	using typename service_t::workItem_t;

public:
	using typename service_t::itemProcessor_t;

	CPUEnergyService();
	CPUEnergyService(CPUEnergyService const&) = delete;
	CPUEnergyService& operator= (CPUEnergyService const&) = delete;
	CPUEnergyService(CPUEnergyService&& copy) {
		*this = std::move(copy);
	}
	CPUEnergyService& operator= (CPUEnergyService&& copy) {
		_d = std::move(copy._d);
		_dataMng = std::move(copy._dataMng);
		return *this;
	}

	virtual ~CPUEnergyService();


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




#endif /* SRC_CPUENERGYSERVICE_H_ */
