/*
 * Worker.tpp
 *
 *  Created on: Mar 26, 2016
 *      Author: uwe
 */

#ifndef SRC_WORKER_TPP_
#define SRC_WORKER_TPP_
#include "/usr/local/cuda/include/cuda_profiler_api.h"
#include "time.h"
//#include "/usr/local/cuda/include/cudaProfiler.h"

namespace as {

template<typename InputType, typename CommonType, typename ResultType>
void Worker<InputType, CommonType, ResultType>::run() {
	bool callAgain = false;
	clock_t start = clock();
	cudaProfilerStart();

	while(true) {
		workItem_t* item;
		if (!callAgain) {
			item = _itemQueue->waitAndPop();
			if (item == nullptr) {
				assert(_itemQueue->terminates());
				break;
			}
		} else {
			bool success = _itemQueue->tryPop(item);
			if (!success) {
				item = nullptr;
			}
		}
		assert(_serviceFnc); // check if fnc object has target
		//
		callAgain = _serviceFnc(item);
		//cudaProfilerStop();
	}
	 cudaProfilerStop();
	 clock_t end = clock();
	 	double elapsed_time = (end - start)/(double)CLOCKS_PER_SEC;
	 	printf("%elapsed time: lf",elapsed_time);
}

} // namespace as

#endif /* SRC_WORKER_TPP_ */
