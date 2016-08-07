/*
 * Int_Service.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef SERVICE_TIMEOUT_H_
#define SERVICE_TIMEOUT_H_

#include <functional>

#include "CPUService.h"

namespace test {

class Service_TimeOut : public as::CPUService<int, int, int> {
public:

	itemProcessor_t createItemProcessor () override;

};

}  // namespace test



#endif /* SERVICE_TIMEOUT_H_ */
