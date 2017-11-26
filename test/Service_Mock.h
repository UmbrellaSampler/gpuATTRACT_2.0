/*
 * Int_Service.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef TEST_TEST_SERVICE_H_
#define TEST_TEST_SERVICE_H_

#include <gmock/gmock.h>

#include "CPUEnergyService.h"
#include "GenericTypes.h"



namespace test {

class Service_Mock : public as::CPUEnergyService<GenericTypes<int, int, int>> {
public:

	MOCK_METHOD0(createItemProcessor, std::function<bool(workItem_t*)> () );

	itemProcessor_t createItemProcessorImpl();

};

}  // namespace test



#endif /* TEST_TEST_SERVICE_H_ */
