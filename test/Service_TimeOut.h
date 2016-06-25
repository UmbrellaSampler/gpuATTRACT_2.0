/*
 * Int_Service.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef SERVICE_TIMEOUT_H_
#define SERVICE_TIMEOUT_H_

#include <functional>

#include "../src/Service.h"

namespace test {

class Service_TimeOut : public as::Service<int, int, int> {
public:

	void initAllocators() override;

	std::function<bool(workItem_t* item)> createItemProcessor () override;

};

}  // namespace test



#endif /* SERVICE_TIMEOUT_H_ */
