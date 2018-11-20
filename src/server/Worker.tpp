/*
 * Worker.tpp
 *
 *  Created on: Mar 26, 2016
 *      Author: uwe
 */

#ifndef SRC_WORKER_TPP_
#define SRC_WORKER_TPP_
#include "cuda_profiler_api.h"
namespace as {

template<typename InputType, typename CommonType, typename ResultType>
void Worker<InputType, CommonType, ResultType>::run() {

	bool callAgain = false;
	//cudaProfilerStart();
	while(true) {
		workItem_t* item;

		if (!callAgain) {
			//printf("size %d id %d\n",  _itemQueue->size(),Thread::id());
			item = _itemQueue->waitAndPop();

			if (item == nullptr) {
			//	printf("terminating %d \n",_itemQueue->terminates());
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
	//	printf("call %d\n",  Thread::id());
		callAgain = _serviceFnc(item);
	}
	//cudaProfilerStop();
}

} // namespace as

#endif /* SRC_WORKER_TPP_ */
