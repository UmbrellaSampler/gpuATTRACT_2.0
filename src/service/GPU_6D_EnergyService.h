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

// ToDo:: make others child of this
template<typename REAL>
class GPU_6D_EnergyService : public GPUService<Types_6D<REAL>> {
public:
	using typename GPUService<Types_6D<REAL>>::service_t;

	using typename TypeWrapper<REAL>::real_t;
	//ToDo: refactor to input_t
	using dof_t = typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;



	/* serves as a builder class*/
//	using configurator_t = Configurator_6D<GPU_6D_EnergyService<REAL>>;
	using configurator_t = Configurator_6D<service_t>;

	/* need public for instantiate the server in configurator */
	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;

//	GPU_6D_EnergyService(std::vector<int> const& deviceIds);
	virtual ~GPU_6D_EnergyService() {};

	distributor_t createDistributor() override;

	itemProcessor_t createItemProcessor() override;

	void setDataManager(std::shared_ptr<DataManager> mng) {
		_dataMng = mng;
	}
private:
	class Private;
	std::shared_ptr<DataManager> _dataMng;

	size_t _workerId; // serves as counter for
	std::vector<int> _deviceIds;


	struct StageResource;

	StageResource createStageResource(workItem_t* item, unsigned const& deviceId);

};

}  // namespace as

#endif


#endif /* SRC_GPU_6D_ENERGYSERVICE_H_ */
