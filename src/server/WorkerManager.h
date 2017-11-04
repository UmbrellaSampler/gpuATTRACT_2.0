/*
 * WorkerManager.h
 *
 *  Created on: Mar 26, 2016
 *      Author: uwe
 */

#ifndef SRC_WORKERMANAGER_H_
#define SRC_WORKERMANAGER_H_

#include <vector>
#include <memory>

#include <Worker.h>
#include "publicTypes.h"

namespace as {

template<typename GenericTypes>
class Service;

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template<typename InputType, typename CommonType, typename ResultType>
class WorkerManager {
	using service_t = Service<GenericTypes<InputType, CommonType, ResultType>>;
	using worker_t = Worker<InputType, CommonType, ResultType>;
	using workItem_t = WorkItem<InputType, CommonType, ResultType>;
public:
	WorkerManager();
	~WorkerManager();

	void setService(service_t* service) noexcept {
		_service = service;
	}

	worker_t* worker(int id) {
		return &_workerPool[id];
	}

	ThreadSafeQueue<workItem_t*>* queue(workerId_t id) noexcept {
		return &_workerQueues[id];
	}

	size_t poolSize() {
		return _workerPool.size();
	}

	size_t queueSize(workerId_t id) {
		return _workerQueues[id].size();
	}


	void pushItemToQueue(workItem_t* item, workerId_t id) {
		_workerQueues[id].push(item);
	}


	void createWorkers(unsigned num);
	void shutDownWorkers() noexcept;

private:
	service_t* _service;
	std::vector<worker_t> _workerPool;
	std::vector<ThreadSafeQueue<workItem_t*>> _workerQueues;

};

}



#endif /* SRC_WORKERMANAGER_H_ */
