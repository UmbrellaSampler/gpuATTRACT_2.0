/*
 * GPU_6D_EnergyService.h
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICE6D_H_
#define SRC_SERVICE_GPUENERGYSERVICE6D_H_

#ifdef CUDA

#include <memory>
#include "GPUEnergyService.h"
#include "publicTypes.h"
#include "Types_6D.h"
#include "nativeTypesWrapper.h"

namespace as {

class DataManager;

template<typename SERVICE>
class Configurator_6D;

template<typename REAL>
class GPUEnergyService6D : public GPUEnergyService<Types_6D<REAL>> {
public:
	using typename GPUEnergyService<Types_6D<REAL>>::service_t;

	using real_t = typename TypeWrapper<REAL>::real_t;
	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;
	using configurator_t = Configurator_6D<REAL>;

	/* need public for instantiate the server in configurator */
	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;


	explicit GPUEnergyService6D(std::shared_ptr<DataManager> dataMng, std::vector<int> const& deviceIds);
	virtual ~GPUEnergyService6D() {};

	distributor_t createDistributor() override;

	itemProcessor_t createItemProcessor() override;

private:
	class Private;

	size_t _workerId; // serves as counter for
	std::vector<int> _deviceIds;


	struct StageResource;

	StageResource createStageResource(workItem_t* item, unsigned const& deviceId);

};

}  // namespace as

#endif


#endif /* SRC_SERVICE_GPUENERGYSERVICE6D_H_ */
