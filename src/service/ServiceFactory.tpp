/*
 * ServiceFactory.tpp
 *
 *  Created on: Aug 12, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_SERVICEFACTORY_TPP_
#define SRC_SERVICE_SERVICEFACTORY_TPP_

#include "ServiceFactory.h"
#include "Service.h"
#include "CmdArgs.h"
#include "CPUEnergyService6D.h"
#include "GenericTypes.h"
#include "CPUEnergyService6DModes.h"
#ifdef CUDA
#include "GPUEnergyService6D.h"
#include "GPUEnergyService6DModes.h"
#include "GPUEnergyServiceMBModes.h"
#endif

namespace as {

//template<typename REAL, template <typename REAL> class GenericTypes>
template<typename REAL>
std::shared_ptr<void> ServiceFactory::create(ServiceType serviceType,
		std::shared_ptr<DataManager> dataMng, CmdArgs const& args)
{
#ifndef CUDA
	(void) args;
#endif
	switch(serviceType) {
		case ServiceType::CPUEnergyService6D:
			return std::shared_ptr<void>( new CPUEnergyService6D<REAL>(dataMng));
		case ServiceType::CPUEnergyService6DModes:
			return std::shared_ptr<void>( new CPUEnergyService6DModes<REAL>(dataMng));
//		case ServiceType::CPUEnergyServiceMBModes:
//			return std::shared_ptr<void>( new CPUEnergyServiceMBModes<REAL>(dataMng));
#ifdef CUDA
		case ServiceType::GPUEnergyService6D:
			return std::shared_ptr<void>( new GPUEnergyService6D<REAL>(dataMng, args.deviceIds));
		case ServiceType::GPUEnergyService6DModes:
			return std::shared_ptr<void>( new GPUEnergyService6DModes<REAL>(dataMng, args.deviceIds));
		case ServiceType::GPUEnergyServiceMBModes:
			return std::shared_ptr<void>( new GPUEnergyServiceMBModes<REAL>(dataMng, args.deviceIds));
#endif
		default:
			throw std::invalid_argument("unknown app to create: " + static_cast<int>(serviceType));
	}
}

} // namespace
#endif /* SRC_SERVICE_SERVICEFACTORY_TPP_ */
