/*
 * ServiceWithDataManager.h
 *
 *  Created on: May 22, 2016
 *      Author: uwe
 */

#ifndef SRC_SERVICEWITHDATAMANAGER_H_
#define SRC_SERVICEWITHDATAMANAGER_H_

#include "Service.h"

namespace as {

/**
 * All services that want to use the data manager should inherit from this class
 */

template<typename InputType, typename CommonType, typename ResultType, typename REAL>
class ServiceWithDataManager : public Service<InputType, CommonType, ResultType> {
	ServiceWithDataManager() {}
	virtual ~ServiceWithDataManager() {}
private:
	// ToDo: Implementation of setters and getters
	std::shared_ptr<DataManager<REAL>> dataMng;
};

}  // namespace as



#endif /* SRC_SERVICEWITHDATAMANAGER_H_ */
