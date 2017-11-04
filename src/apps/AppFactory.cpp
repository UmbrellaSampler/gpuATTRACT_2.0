/*
 * AppFactory.cpp
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include <exception>
#include "AppFactory.h"
#include "scATTRACT.h"
#include "CmdArgs.h"
#include <CPUEnergyService6D.h>

#ifdef CUDA
#include <GPUEnergyService6D.h>
#endif

namespace as {

std::unique_ptr<App> AppFactory::create(const CmdArgs& args) {

	std::unique_ptr<App> app;

	Platform p = Platform::unspecified;
	if (args.numCPUs > 0) {
		p = Platform::CPU;
	} else if (args.deviceIds.size() > 0) {
		p = Platform::GPU;
	}

	if (args.precision == "single") {
		app = create<float>(args.app, p);
	} else if (args.precision == "double") {
		app = create<double>(args.app, p);
	} else {
		throw std::invalid_argument("unknown precision specification: " + args.precision);
	}

	return app;
}

template<typename REAL>
std::unique_ptr<App> AppFactory::create(AppType appType, Platform p) {
//std::unique_ptr<App> AppFactory::create(AppType appType, ServiceType ServiceType, Platform p) {
// for HM and/or MultiBodies --> need TypeFactory (like TypeWrapper) for EnergySerives

	std::unique_ptr<App> app;

	switch (appType) {
	case AppType::Score:

		switch (p) {
			case Platform::CPU:
				app = std::unique_ptr<App> (new scATTRACT<CPUEnergyService6D<REAL>>());
				break;
#ifdef CUDA
			case Platform::GPU:
				app = std::unique_ptr<App> (new scATTRACT<GPUEnergyService6D<REAL>>());
				break;
#endif
			default:
				throw std::invalid_argument("unknown platform to create: " + static_cast<int>(appType));
		}

		// in case of multiple possible Service types (GPU, Modes) create private helper template factory method
		break;
	default:
		throw std::invalid_argument("unknown AppType: " + static_cast<int>(appType));
	}

	return app;
}

} // namespace as

