/*
 * Worker.h
 *
 *  Created on: Mar 26, 2016
 *      Author: uwe
 */

#ifndef SRC_WORKER_H_
#define SRC_WORKER_H_

#include "ThreadSafeQueue.h"
//#include "cuda_profiler_api.h"
namespace as {

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template<typename InputType, typename CommonType, typename ResultType>
class Worker : public Thread {
	using workItem_t = WorkItem<InputType, CommonType, ResultType>;
public:

	void setItemQueue(ThreadSafeQueue<workItem_t*>* queue ) noexcept {
		_itemQueue = queue;
	}

	void setFncObj(const std::function<bool(workItem_t*)>& obj) {
		_serviceFnc = obj;
	}

	void signalTerminate () {
		_itemQueue->signalTerminate();
	}

private:

	void run() override;

	std::function<bool(workItem_t*)> _serviceFnc;
	ThreadSafeQueue<workItem_t*>* _itemQueue;

};

} // namespace



#endif /* SRC_WORKER_H_ */
