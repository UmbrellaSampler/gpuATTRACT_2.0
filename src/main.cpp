/*
 * main.cpp
 *
 *  Created on: Dec 29, 2015
 *      Author: uwe
 */

#include <iostream>
#include <numeric>

#include "ServerIncludes.h"
#include "Server.h"
#include "Request.h"
#include "Allocator.h"
#include "TestService.h"
#include "macros.h"


using namespace std;

int main(int argc, char* argv[]) {
	as::Server<TestService> server;
	std::vector<unsigned> itemSizes = {1,2,3,5,7,55,111,390,557,1000};
	std::vector<unsigned> inputSizes = {10,537,1001,1,2,3,2000, 3000, 4000, 100000};
	std::vector<unsigned> poolSize = {1,2,3,4};

	for (unsigned workers : poolSize) {
		server.createWorkers(workers);
		cout << workers << " Workers" << endl;

		for (unsigned inputSize : inputSizes) {
			int input[inputSize];
			int result[inputSize];
			std::iota(input, input + inputSize, 0);
			float common = 7.459f;
			as::Request<int, float> request;
			for (unsigned itemSize : itemSizes) {
				request.setDataAndSize(input, inputSize);
				request.setCommon(common);
				server.setItemSize(itemSize);
				server.submit(request);
				server.wait(request, result);
				for(unsigned i = 0; i < inputSize; ++i) {
					ASSERT(result[i] == static_cast<int>(input[i]*common));
				}
			}
		}
	}
	cout << "exit=OK" << endl;

	return 0;

}

