/*
 * Int_Service.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */
#include <memory>

#include "Service_Mock.h"
#include "../src/Allocator.h"

#include "../src/WorkItem.h"

std::function<bool(test::Service_Mock::workItem_t* item)> test::Service_Mock::createItemProcessor_fake() {

	std::function<bool(workItem_t* item)> fncObj = [] (workItem_t* item) {
		using namespace std;
		auto size = item->size();
		auto input = item->inputBuffer();
		auto result = item->resultBuffer();
		auto common = *item->common();
		for(unsigned i = 0; i < size; ++i) {
			result[i] = input[i]*common;
		}
		for(unsigned i = 0; i < size; ++i) {
			for (unsigned j = 0; j < 13; ++j) {
				float div = common/2.0f;
				result[i] /= div;
			}
		}
		item->setProcessed();
		return false;
	};

	return fncObj;
}

void test::Service_Mock::initAllocators() {
	setInputAllocator(std::make_shared<as::HostAllocator<input_t>>());
	setResultAllocator(std::make_shared<as::HostAllocator<result_t>>());
}


