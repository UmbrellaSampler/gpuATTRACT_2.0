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

#ifdef CUDA
#include "GPUEnergyService6D.h"
#endif

namespace as {

template<typename REAL, template <typename REAL> class GenericTypes>
std::unique_ptr<Service<GenericTypes<REAL>>> ServiceFactory::create(ServiceType serviceType,
		std::shared_ptr<DataManager> dataMng, CmdArgs const& args)
{
#ifndef CUDA
	(void) args;
#endif
	switch(serviceType) {
		case ServiceType::CPUEnergyService6D:
			return std::unique_ptr<Service<GenericTypes<REAL>>>( new CPUEnergyService6D<REAL>(dataMng));
#ifdef CUDA
		case ServiceType::GPUEnergyService6D:
			return std::unique_ptr<Service<GenericTypes<REAL>>>( new GPUEnergyService6D<REAL>(dataMng, args.deviceIds));
#endif
		default:
			throw std::invalid_argument("unknown app to create: " + static_cast<int>(serviceType));
	}
}

} // namespace
#endif /* SRC_SERVICE_SERVICEFACTORY_TPP_ */
