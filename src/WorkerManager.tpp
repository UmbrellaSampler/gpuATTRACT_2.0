/*
 * WorkerManager.tpp
 *
 *  Created on: Mar 26, 2016
 *      Author: uwe
 */

#ifndef SRC_WORKERMANAGER_TPP_
#define SRC_WORKERMANAGER_TPP_

#include <queue>
#include "WorkerManager.h"
#include "Service.h"


using namespace as;

template<typename InputType, typename CommonType, typename ResultType>
WorkerManager<InputType, CommonType, ResultType>::WorkerManager() :
	_service(nullptr)
{}

template<typename InputType, typename CommonType, typename ResultType>
WorkerManager<InputType, CommonType, ResultType>::~WorkerManager() {
	shutDownWorkers();
}

template<typename InputType, typename CommonType, typename ResultType>
void WorkerManager<InputType, CommonType, ResultType>::shutDownWorkers() noexcept{
	for (auto& worker : _workerPool) {
		worker.signalTerminate();
		worker.join();
	}
	for (auto const& queue : _workerQueues) {
		assert(queue.empty());
	}
}

template<typename InputType, typename CommonType, typename ResultType>
void WorkerManager<InputType, CommonType, ResultType>::createWorkers(unsigned num) {
	assert(_service != nullptr);
	if (_workerPool.size() > 0) {
		shutDownWorkers();
	}
	_workerPool = std::move(std::vector<worker_t>(num));
	_workerQueues = std::move(std::vector<ThreadSafeQueue<workItem_t*>>(num));
	for (int i = 0; i < num; ++i) {
		auto& worker = _workerPool[i];
		worker.setFncObj(_service->createItemProcessor());
		auto& queue = _workerQueues[i];
		worker.setItemQueue(&queue);
		worker.start();
	}
}



#endif /* SRC_WORKERMANAGER_TPP_ */
