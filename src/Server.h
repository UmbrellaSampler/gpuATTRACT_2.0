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

template<class T>
class ThreadSafeQueue;

constexpr unsigned defaultItemSize = 100;
constexpr unsigned maxWaitMilliSeconds = 20000;

template<typename Service>
class Server {
public:
	using input_t = typename Service::input_t;
	using result_t = typename Service::result_t;
	using common_t = typename Service::common_t;

	using workItem_t = typename Service::workItem_t;

	using request_t = Request<input_t, common_t>;

	Server();
	explicit Server(std::shared_ptr<Service> service);
	Server(Server const&) = default;

	~Server();

	std::shared_ptr<Service> service() noexcept {
		return _service;
	}

	void setService(std::shared_ptr<Service> const& service);

	void createWorkers(unsigned number);

	void setWaitTime(unsigned milliSeconds) noexcept {
		_waitTime = milliSeconds;
	}

	void setItemSize(size_t itemSize);

	size_t itemSize() const noexcept{
		return _itemSize;
	}

	void submit(request_t& req);

	void wait(request_t const& req, result_t* result);

private:

	void configureBufferAllocators();
	void attachServerBuffers(request_t const*);
	void copyRequestBuffer(request_t const*);
	void createWorkItemsAndPush(request_t const*);
	void synchronizeWith(request_t const*);
	void returnServerBuffers(request_t const*);
	void copyResultBuffer(request_t const*, result_t*);

	std::shared_ptr<Service> _service;
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
