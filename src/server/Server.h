/*
 * Server.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef SERVER_H_
#define SERVER_H_

#include <memory>
#include <set>

#include "Service.h"

namespace as {

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template<typename InputType, typename CommonType>
class Request;

template<typename InputType, typename ResultType>
class BufferManager;

template<typename InputType, typename CommonType, typename ResultType>
class RequestManager;

template<typename InputType, typename CommonType, typename ResultType>
class Dispatcher;

template<typename InputType, typename CommonType, typename ResultType>
class WorkerManager;

template<typename T>
class ThreadSafeQueue;

constexpr unsigned DEFAULT_ITEM_SIZE = 100;
constexpr unsigned MAX_WAIT_MILLISECONDS = 20000;

template<typename GenericTypes>
class Server {
public:
	using service_t = Service<GenericTypes>;
	using input_t = typename GenericTypes::input_t;
	using result_t = typename GenericTypes::result_t;
	using common_t = typename GenericTypes::common_t;

	using workItem_t = typename service_t::workItem_t;

	using request_t = Request<input_t, common_t>;

//	Server();
	explicit Server(std::shared_ptr<service_t> service);
	Server(Server const&) = default;

	~Server();

	std::shared_ptr<service_t> service() noexcept {
		return _service;
	}

	void setService(std::shared_ptr<service_t> const& service);

	void createWorkers(unsigned number);

	void setWaitTime(unsigned milliSeconds) noexcept {
		_waitTime = milliSeconds;
	}

	void setItemSize(size_t itemSize);

	size_t itemSize() const noexcept{
		return _itemSize;
	}

	void submit(std::shared_ptr<request_t> const& request);

	void wait(std::shared_ptr<request_t> const& request, result_t* result);

private:

	void configureBufferAllocators();
	void attachServerBuffers(std::shared_ptr<request_t> const&);
	void copyRequestBuffer(std::shared_ptr<request_t> const&);
	void createWorkItemsAndPush(std::shared_ptr<request_t> const&);
	void synchronizeWith(std::shared_ptr<request_t> const& request);
	void returnServerBuffers(std::shared_ptr<request_t> const&);
	void copyResultBuffer(std::shared_ptr<request_t> const&, result_t*);

	std::shared_ptr<service_t> _service;
	size_t _itemSize;
	unsigned _waitTime;

	std::unique_ptr<BufferManager<input_t, result_t>> _bufferMng;
	std::unique_ptr<RequestManager<input_t, common_t, result_t>> _requestMng;
	std::unique_ptr<Dispatcher<input_t, common_t, result_t>> _dispatcher;
	std::unique_ptr<WorkerManager<input_t, common_t, result_t>> _workerMng;
	std::unique_ptr<ThreadSafeQueue<std::vector<workItem_t>*>> _requestQueue;


};

}  // namespace as



#endif /* SERVER_H_ */
