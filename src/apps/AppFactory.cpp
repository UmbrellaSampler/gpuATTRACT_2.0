/*
 * AppFactory.cpp
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include <exception>
#include "AppFactory.h"
#include "scATTRACT.h"
#include "CPU_6D_EnergyService.h"

#ifdef CUDA
#include "GPU_6D_EnergyService.h"
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

	std::unique_ptr<App> app;

	switch (appType) {
	case AppType::Score:

		switch (p) {
			case Platform::CPU:
				app = std::unique_ptr<App> (new scATTRACT<CPU_6D_EnergyService<REAL>>());
				break;
#ifdef CUDA
			case Platform::GPU:
				app = std::unique_ptr<App> (new scATTRACT<GPU_6D_EnergyService<REAL>>());
				break;
#endif
			default:
				throw std::invalid_argument("unknown platform to create: " + static_cast<int>(appType));
		}

		// in case of multiple possible Service types (GPU, Modes) create private helper template factory method
		break;
	default:
		throw std::invalid_argument("unknown app to create: " + static_cast<int>(appType));
	}

	return app;
}

} // namespace as

