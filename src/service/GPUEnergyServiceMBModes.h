/*
 * GPU_MB_EnergyService.h
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICE_GPUENERGYSERVICEMBMODES_H_
#define SRC_SERVICE_GPUENERGYSERVICEMBMODES_H_

#ifdef CUDA

#include <memory>
#include "GPUEnergyService.h"
#include "publicTypes.h"
#include "Types_MB_Modes.h"
#include "nativeTypesWrapper.h"

namespace as {

class DataManager;

template<typename SERVICE>
class Configurator_MB_Modes;

template<typename REAL>
class GPUEnergyServiceMBModes : public GPUEnergyService<Types_MB_Modes<REAL>> {
public:
	using typename GPUEnergyService<Types_MB_Modes<REAL>>::service_t;

	using real_t = typename TypeWrapper<REAL>::real_t;
	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;
	using configurator_t = Configurator_MB_Modes<REAL>;

	/* need public for instantiate the server in configurator */
	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;
	using typename service_t::distributor_t;


	explicit GPUEnergyServiceMBModes(std::shared_ptr<DataManager> dataMng, std::vector<int> const& deviceIds);
	virtual ~GPUEnergyServiceMBModes() {};

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


#endif /* SRC_SERVICE_GPUENERGYSERVICEMB_H_ */
