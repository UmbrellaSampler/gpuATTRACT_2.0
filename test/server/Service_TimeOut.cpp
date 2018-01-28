/*
 * Int_Service.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */
#include <memory>
#include "Service_TimeOut.h"
#include "Allocator.h"

#include "WorkItem.h"

auto test::Service_TimeOut::createItemProcessor() -> itemProcessor_t{

	std::function<bool(workItem_t* item)> fncObj = [] (workItem_t* item) {
		return false;
	};

	return fncObj;
}




