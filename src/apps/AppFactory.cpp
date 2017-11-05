/*
 * AppFactory.cpp
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include <exception>
#include "AppFactory.h"
#include "scATTRACT.h"
#include "emATTRACT.h"
#include "CmdArgs.h"

namespace as {

std::unique_ptr<App> AppFactory::create(const CmdArgs& args) {

	std::unique_ptr<App> app;
	if (args.precision == "single") {
		app = create<float>(args.app);
	} else if (args.precision == "double") {
		app = create<double>(args.app);
	} else {
		throw std::invalid_argument("unknown precision specification: " + args.precision);
	}

	return app;
}

template<typename REAL>
std::unique_ptr<App> AppFactory::create(AppType appType) {
//std::unique_ptr<App> AppFactory::create(AppType appType, ServiceType ServiceType, Platform p) {
// for HM and/or MultiBodies --> need TypeFactory (like TypeWrapper) for EnergySerives

	std::unique_ptr<App> app;

	switch (appType) {
	case AppType::SCORE:
		app = std::unique_ptr<App> (new scATTRACT<Types_6D<REAL>>());
		break;

	case AppType::EM:
		app = std::unique_ptr<App> (new emATTRACT<Types_6D<REAL>>());
		break;

	default:
		throw std::invalid_argument("unknown AppType: " + static_cast<int>(appType));
	}

	return app;
}

} // namespace as

