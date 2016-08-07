/*
 * Dispatcher.h
 *
 *  Created on: Mar 24, 2016
 *      Author: uwe
 */

#ifndef SRC_DISPATCHER_H_
#define SRC_DISPATCHER_H_

#include "Service.h"
#include "Thread.h"
#include "publicTypes.h"

namespace as {

template<typename InputType, typename CommonType, typename ResultType>
class RequestManager;

template<typename InputType, typename CommonType, typename ResultType>
class WorkerManager;

template<typename InputType, typename CommonType>
class Request;

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template <class T>
class ThreadSafeQueue;

struct FillLevel;

template<typename InputType, typename CommonType, typename ResultType>
class Dispatcher : public Thread {
	using workItem_t = WorkItem<InputType, CommonType, ResultType>;
	using requestManager_t = RequestManager<InputType, CommonType, ResultType>;
	using workerManager_t = WorkerManager<InputType, CommonType, ResultType>;
	using distributor_t = typename Service<InputType, CommonType, ResultType>::distributor_t;
public:

	Dispatcher() :
		_workerMng(nullptr),
		_requestQueue(nullptr),
		_getWorkerIds(nullptr)
	{};

	~Dispatcher() {};

	void setWorkerManager(workerManager_t* workerMng) noexcept {
		_workerMng = workerMng;
	}

	void setQueue(ThreadSafeQueue<std::vector<workItem_t>*>* queue) noexcept {
		_requestQueue = queue;
	}

	void setDistributor(std::function<std::vector<as::workerId_t>(CommonType const*, size_t numWorkers)> distributor) {
		_getWorkerIds = distributor;
	}

	void signalTerminate() noexcept {
		_requestQueue->signalTerminate();
	}


private:
	void run() override;
	void dispatch(std::vector<workItem_t>* items);
	std::set<unsigned> getWorkerIds(CommonType const* common) const noexcept;
	std::vector<FillLevel> getFillLevelsSorted(const std::vector<unsigned>& workers) const noexcept;
	void distributeItem(workItem_t&, std::vector<FillLevel>&);

	workerManager_t* _workerMng;
	ThreadSafeQueue<std::vector<workItem_t>*>* _requestQueue;

	distributor_t _getWorkerIds;

};

}  // namespace as



#endif /* SRC_DISPATCHER_H_ */
