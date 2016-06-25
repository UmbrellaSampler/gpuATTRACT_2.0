/*
 * RequestManager.cpp
 *
 *  Created on: Mar 21, 2016
 *      Author: uwe
 */
#ifndef SRC_REQUESTMANAGER_TPP_
#define SRC_REQUESTMANAGER_TPP_

#include <cassert>
#include <atomic>
#include <mutex>

#include "RequestManager.h"
#include "Request.h"
#include "ServerBuffers.h"
#include "RequestMembers.h"


using namespace as;

template<typename InputType, typename CommonType, typename ResultType>
void RequestManager<InputType, CommonType, ResultType>::registerRequest(request_t const* request) {
	assert(!isValid(request));
	auto& requestMembers = _requestContainer[request]; // implicit construction of new entry

	if (request->size() > 0) {
		auto serverBuffers = requestMembers.serverBuffers();
		serverBuffers->size = request->size();

		CommonType* common = requestMembers.common();
		*common = request->common();
	}
}

template<typename InputType, typename CommonType, typename ResultType>
void RequestManager<InputType, CommonType, ResultType>::removeRequest(request_t const* request) {
	_requestContainer.erase(request);
}

template<typename InputType, typename CommonType, typename ResultType>
bool RequestManager<InputType, CommonType, ResultType>::isValid(request_t const* request) const {
	return _requestContainer.find(request) != _requestContainer.end();
}

template<typename InputType, typename CommonType, typename ResultType>
size_t RequestManager<InputType, CommonType, ResultType>::containerSize() const {
	return _requestContainer.size();
}

template<typename InputType, typename CommonType, typename ResultType>
ServerBuffers<InputType, ResultType>* RequestManager<InputType, CommonType, ResultType>::buffers(request_t const* request) {
	auto& requestMembers = _requestContainer.at(request);
	return requestMembers.serverBuffers();
}

template<typename InputType, typename CommonType, typename ResultType>
std::vector<WorkItem<InputType, CommonType, ResultType>>* RequestManager<InputType, CommonType, ResultType>::workItems(request_t const* request) {
	auto& requestMembers = _requestContainer.at(request);
	return requestMembers.items();
}

template<typename InputType, typename CommonType, typename ResultType>
bool RequestManager<InputType, CommonType, ResultType>::isProcessed(request_t const* request) {
	auto items = workItems(request);
	for (auto& item : *items) {
		if (!item.isProcessed()) {
			return false;
		}
	}
	std::atomic_thread_fence(std::memory_order_acquire);
	return true;
}

template<typename InputType, typename CommonType, typename ResultType>
CommonType* RequestManager<InputType, CommonType, ResultType>::common(request_t const* request) {
	auto& requestMembers = _requestContainer.at(request);
	return requestMembers.common();
}

#endif /* SRC_REQUESTMANAGER_TPP_ */
