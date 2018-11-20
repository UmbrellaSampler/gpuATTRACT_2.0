/*
 * ServiceFactory.h
 *
 *  Created on: Aug 12, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_SERVICEFACTORY_H_
#define SRC_SERVICE_SERVICEFACTORY_H_

#include <memory>
#include "ServiceType.h"
#include "GenericTypes.h"

namespace as {

template<typename GenericTypes>
class Service;

class CmdArgs;
class DataManager;

class ServiceFactory {
public:
	ServiceFactory() = delete;

//	template<typename REAL, template <typename REAL> class GenericTypes>
//	static std::unique_ptr<Service<GenericTypes<REAL>>> create(ServiceType, std::shared_ptr<DataManager>, CmdArgs const&);

	template<typename REAL>
	static std::shared_ptr<void> create(ServiceType, std::shared_ptr<DataManager>, CmdArgs const&, int threadsPerDevice);

};

}  // namespace as

#endif /* SRC_SERVICE_SERVICEFACTORY_H_ */
