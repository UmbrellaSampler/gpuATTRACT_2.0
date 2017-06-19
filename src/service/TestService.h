/*
 * TestService.h
 *
 *  Created on: Mar 22, 2016
 *      Author: uwe
 */

#ifndef SRC_TESTSERVICE_H_
#define SRC_TESTSERVICE_H_

#include "Service.h"
#include <functional>
#include <cmath>

class TestService : public as::Service<int, float, int> {
public:

	std::function<bool(workItem_t* item)> createItemProcessor() {
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
				for (unsigned j = 0; j < 1000; ++j) {
					float div = 1.0;
					div *= common;
					div /= common;
					result[i] /= std::round(div);
				}
			}
			item->setProcessed();
			return false;
		};

		return fncObj;
	}

	void initAllocators() {
		setInputAllocator(std::make_shared<as::HostAllocator<input_t>>());
		setResultAllocator(std::make_shared<as::HostAllocator<result_t>>());
	}
};

#endif /* SRC_TESTSERVICE_H_ */
