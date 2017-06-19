/*
 * GPU_6D_EnergyService.h
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_GPU_6D_ENERGYSERVICE_H_
#define SRC_GPU_6D_ENERGYSERVICE_H_

#ifdef CUDA

#include <memory>
#include "GPUService.h"
#include "publicTypes.h"
#include "Types_6D.h"
#include "nativeTypesWrapper.h"

namespace as {

class DataManager;

template<typename SERVICE>
class Configurator_6D;

template<typename REAL>
class GPU_6D_EnergyService : public GPUService<typename Types_6D<REAL>::DOF, typename Types_6D<REAL>::Common,  typename Types_6D<REAL>::Result> {
public:
	using real_t = typename TypeWrapper<REAL>::real_t;
	using dof_t = typename Types_6D<real_t>::DOF;
	using common_t = typename Types_6D<real_t>::Common;
	using result_t = typename Types_6D<real_t>::Result;
public:
	using typename GPUService<dof_t, common_t,  result_t>::distributor_t;
	using service_t = GPUService<dof_t, common_t, result_t>;
public:

	/* serves as a builder class*/
	using configurator_t = Configurator_6D<GPU_6D_EnergyService<REAL>>;

	/* need public for instantiate the server in configurator */

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;

	distributor_t createDistributor() override;

	itemProcessor_t createItemProcessor() override;

	void setDataManager(std::shared_ptr<DataManager> mng) {
		_dataMng = mng;
	}
private:
	class Private;
	std::shared_ptr<DataManager> _dataMng;

	size_t _workerId; // serves as counter for


	struct StageResource;

	StageResource createStageResource(workItem_t* item, unsigned const& deviceId);

};

}  // namespace as

#endif


#endif /* SRC_GPU_6D_ENERGYSERVICE_H_ */
