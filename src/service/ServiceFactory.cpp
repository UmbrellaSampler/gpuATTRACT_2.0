/*
 * ServiceFactory.cpp
 *
 *  Created on: Aug 12, 2017
 *      Author: uwe
 */


#include "ServiceFactory.tpp"
#include "Types_6D.h"
#include "Types_6D_Modes.h"


namespace as {

template
std::shared_ptr<void> ServiceFactory::create<float>(ServiceType serviceType,
		std::shared_ptr<DataManager> dataMng, CmdArgs const& args, int threadsPerDevice =1);

template
std::shared_ptr<void> ServiceFactory::create<double>(ServiceType serviceType,
		std::shared_ptr<DataManager> dataMng, CmdArgs const& args, int threadsPerDevice =1);

//template
//std::unique_ptr<Service<Types_6D_Modes<float>>> ServiceFactory::create<float, Types_6D_Modes>(ServiceType serviceType,
//		std::shared_ptr<DataManager> dataMng, CmdArgs const& args);
//
//template
//std::unique_ptr<Service<Types_6D_Modes<double>>> ServiceFactory::create<double, Types_6D_Modes>(ServiceType serviceType,
//		std::shared_ptr<DataManager> dataMng, CmdArgs const& args);

} // namespace


