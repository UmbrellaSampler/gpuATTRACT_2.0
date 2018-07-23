/*
 * AppFactory.cpp
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include <stdexcept>
#include <type_traits>
#include <iostream>
#include "AppFactory.h"
#include "scATTRACT.h"
#include "emATTRACT.h"
#include "mcATTRACT.h"
#include "CmdArgs.h"


namespace as {

std::unique_ptr<App> AppFactory::create(const CmdArgs& args) {

	std::unique_ptr<App> app;
	const bool useModes = args.numModes > 0;
	if (args.precision == "single") {
		if (useModes) {
			Common_Modes::numModesRec = args.numModes;
			Common_Modes::numModesLig = args.numModes;
			app = create<Types_6D_Modes<float>>(args.app);
		} else {
			app = create<Types_6D<float>>(args.app);
		}
	} else if (args.precision == "double") {
		if (useModes) {
			Common_Modes::numModesRec = args.numModes;
			Common_Modes::numModesLig = args.numModes;
			app = create<Types_6D_Modes<double>>(args.app);

		} else {
			app = create<Types_6D<double>>(args.app);
		}
	} else {
		throw std::invalid_argument("unknown precision specification: " + args.precision);
	}


	return app;
}

template<typename GenericTypes>
std::unique_ptr<App> AppFactory::create(AppType appType) {
//std::unique_ptr<App> AppFactory::create(AppType appType, ServiceType ServiceType, Platform p) {
// for HM and/or MultiBodies --> need TypeFactory (like TypeWrapper) for EnergySerives

	std::unique_ptr<App> app;

	switch (appType) {
	case AppType::SCORE:
		app = std::unique_ptr<App> (new scATTRACT<GenericTypes>());
		break;

	case AppType::EM:
		app = std::unique_ptr<App> (new emATTRACT<GenericTypes>());
		break;
	case AppType::MC:
		app = std::unique_ptr<App> (new mcATTRACT<GenericTypes>());
		break;
 	default:
		throw std::invalid_argument("unknown AppType: " + static_cast<int>(appType));
	}

	return app;
}

} // namespace as

