/*
 * Int_Service.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef SERVICE_TIMEOUT_H_
#define SERVICE_TIMEOUT_H_

#include <functional>

#include "../src/service/CPUEnergyService.h"

namespace test {

class Service_TimeOut : public as::CPUEnergyService<int, int, int> {
public:

	itemProcessor_t createItemProcessor () override;

};

}  // namespace test



#endif /* SERVICE_TIMEOUT_H_ */
