/*
 * ServiceFactory.cpp
 *
 *  Created on: Aug 12, 2017
 *      Author: uwe
 */


#include "ServiceFactory.tpp"
#include "Types_6D.h"

namespace as {

template
std::unique_ptr<Service<Types_6D<float>>> ServiceFactory::create<float, Types_6D>(ServiceType serviceType,
		std::shared_ptr<DataManager> dataMng, CmdArgs const& args);

template
std::unique_ptr<Service<Types_6D<double>>> ServiceFactory::create<double, Types_6D>(ServiceType serviceType,
		std::shared_ptr<DataManager> dataMng, CmdArgs const& args);


} // namespace


