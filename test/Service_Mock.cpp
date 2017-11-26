/*
 * Int_Service.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */
#include <memory>

#include "Service_Mock.h"

#include "WorkItem.h"

auto test::Service_Mock::createItemProcessorImpl() -> itemProcessor_t {

	itemProcessor_t fncObj = [] (workItem_t* item) {
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


