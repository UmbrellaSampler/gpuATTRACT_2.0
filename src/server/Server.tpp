/*
 * Server.cpp
 *
 *  Created on: Dec 29, 2015
 *      Author: uwe
 */


/* this file does not necessarily need to be compiled */

#ifndef SRC_SERVER_TPP_
#define SRC_SERVER_TPP_

#include <vector>
#include <stdexcept>
#include <thread>
#include <chrono>

#include "Server.h"

#include "../commons/helper.h"
#include "Service.h"
#include "WorkItem.h"
#include "Request.h"
#include "BufferManager.h"
#include "RequestManager.h"
#include "Dispatcher.h"
#include "WorkerManager.h"
#include "ThreadSafeQueue.h"
#include "ServerFunctions.h"
#include "macros.h"

namespace as {

template<typename Service>
Server<Service>::Server(std::shared_ptr<Service> service):
	_itemSize(defaultItemSize),
	_waitTime(maxWaitMilliSeconds),
	_bufferMng(make_unique<BufferManager<input_t, result_t>>()),
	_requestMng(make_unique<RequestManager<input_t, common_t, result_t>>()),
	_dispatcher(make_unique<Dispatcher<input_t, common_t, result_t>>()),
	_workerMng(make_unique<WorkerManager<input_t, common_t, result_t>>()),
	_requestQueue(make_unique<ThreadSafeQueue<std::vector<workItem_t>*>>())
{
	try {
		setService(service);
	} catch (std::exception& e) {
		throw;
	}

	_dispatcher->setWorkerManager(_workerMng.get());
	_dispatcher->setQueue(_requestQueue.get());
	_dispatcher->start();
};

template<typename Service>
Server<Service>::Server():
	Server(([]() -> std::shared_ptr<Service> {
			std::shared_ptr<Service> service = std::make_shared<Service>();
			service->initAllocators();
			return service;
		})()
	)
{};

template<typename Service> // needed for invalid application of 'sizeof' compiler error messages
Server<Service>::~Server()
{
	_dispatcher->signalTerminate();
	_dispatcher->join();

};

template<typename Service>
void Server<Service>::setService(std::shared_ptr<Service> const& service) {
	if (service.get() == nullptr) {
		throw std::invalid_argument("Invalid service (nullptr).");
	}
	_service = service;
	try {
		configureBufferAllocators();
	} catch (std::exception& e) {
		throw;
	}
	_workerMng->setService(_service.get());
	_dispatcher->setDistributor(_service->createDistributor());

}

template<typename Service>
void Server<Service>::createWorkers(unsigned number) {
	assert(_service != nullptr);
	_workerMng->createWorkers(number);
}

template<typename Service>
void Server<Service>::setItemSize(size_t itemSize) {
	if(itemSize == 0) {
		throw std::invalid_argument("Item size must be greater than zero.");
	}
	_itemSize = itemSize;
}

template<typename Service>
void Server<Service>::submit(request_t& request) {

	/* Client */
	{
		if (_requestMng->isValid(&request)) {
			throw std::invalid_argument("Request has already been submitted.");
		} else if (request.data() == nullptr) {
			throw std::invalid_argument("Invalid input buffer (nullptr).");
		}

		_requestMng->registerRequest(&request);

		if(request.size() > 0) {
			attachServerBuffers(&request);
			copyRequestBuffer(&request);
			createWorkItemsAndPush(&request);
		}

		request.reset();
	}

}

template<typename Service>
void Server<Service>::attachServerBuffers(request_t const* request) {
	ServerBuffers<input_t, result_t>* buffers = _requestMng->buffers(request);
	auto& bufferSize = buffers->size;
	buffers->inputBuffer = _bufferMng->getInputBuffer(bufferSize);
	buffers->resultBuffer = _bufferMng->getResultBuffer(bufferSize);
}

template<typename Service>
void Server<Service>::copyRequestBuffer(request_t const* request) {
	ServerBuffers<input_t, result_t>* buffers = _requestMng->buffers(request);
	const auto& bufferSize = buffers->size;
	assert(bufferSize == request->size());
	input_t *const origBuffer = request->data();
	input_t *const serverBuffer = buffers->inputBuffer;
	std::copy(origBuffer, origBuffer + bufferSize, serverBuffer);
}

template<typename Service>
void Server<Service>::createWorkItemsAndPush(request_t const* request) {
	ServerBuffers<input_t, result_t> const* buffers = _requestMng->buffers(request);
	auto const& requestSize = buffers->size;
	common_t* common = _requestMng->common(request);
	std::vector<workItem_t>* items = _requestMng->workItems(request);
	*items = std::move(createWorkItemsFn(buffers, requestSize, common, _itemSize));
	_requestQueue->push(items);
}



template<typename Service>
void Server<Service>::wait(request_t const& request, result_t* clientBuffer)
{
	if (!_requestMng->isValid(&request) || clientBuffer == nullptr) {
		throw std::invalid_argument("Either request has not been submitted or clientBuffer invalid (nullptr).");
	}

	try {
		synchronizeWith(&request);
	} catch (std::exception& e) {
		returnServerBuffers(&request);
		_requestMng->removeRequest(&request);
		throw;
	}

	copyResultBuffer(&request, clientBuffer);
	returnServerBuffers(&request);
	_requestMng->removeRequest(&request);
}

template<typename Service>
void Server<Service>::synchronizeWith(request_t const* request) {
	unsigned count = 0;
	bool processed;
	while (!(processed = _requestMng->isProcessed(request)) && count < _waitTime) {
		++count;
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}

	if (!processed) {
		throw std::runtime_error("Request takes too long to finish for unknown reasons");
	}
}

template<typename Service>
void Server<Service>::returnServerBuffers(request_t const* request) {
	ServerBuffers<input_t, result_t> const* buffers = _requestMng->buffers(request);
	auto const& size = buffers->size;
	if (size > 0) {
		_bufferMng->returnInputBuffer(buffers->inputBuffer);
		_bufferMng->returnResultBuffer(buffers->resultBuffer);
	}
}

template<typename Service>
void Server<Service>::copyResultBuffer(request_t const* request, result_t* clientBuffer) {
	ServerBuffers<input_t, result_t> const* buffers = _requestMng->buffers(request);
	auto const& size = buffers->size;
	if (size > 0) {
		const result_t* resultBuffer = buffers->resultBuffer;
		std::copy(resultBuffer, resultBuffer + size, clientBuffer);
	}
}

template<typename Service>
void Server<Service>::configureBufferAllocators() {
	if (_service->getInputAllocator() == nullptr ||
		_service->getResultAllocator() == nullptr) {
		throw std::invalid_argument("Invalid service allocators (nullptr).");
	}
	_bufferMng->setInputAllocator(_service->getInputAllocator());
	_bufferMng->setResultAllocator(_service->getResultAllocator());
}

} // namespace as

#endif /* SRC_SERVER_TPP_ */

