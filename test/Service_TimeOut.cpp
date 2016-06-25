/*
 * Int_Service.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */
#include <memory>
#include "Service_TimeOut.h"
#include "../src/Allocator.h"

#include "../src/WorkItem.h"

std::function<bool(test::Service_TimeOut::workItem_t* item)> test::Service_TimeOut::createItemProcessor() {

	std::function<bool(workItem_t* item)> fncObj = [] (workItem_t* item) {
		return false;
	};

	return fncObj;
}

void test::Service_TimeOut::initAllocators() {
	setInputAllocator(std::make_shared<as::HostAllocator<input_t>>());
	setResultAllocator(std::make_shared<as::HostAllocator<result_t>>());
}


