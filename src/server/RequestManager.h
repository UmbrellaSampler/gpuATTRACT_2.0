/*
 * RequestManager.h
 *
 *  Created on: Mar 21, 2016
 *      Author: uwe
 */

#ifndef SRC_REQUESTMANAGER_H_
#define SRC_REQUESTMANAGER_H_

#include <set>
#include <map>
#include "ThreadSafeQueue.h"
#include <memory>
namespace as {

template<typename InputType, typename CommonType>
class Request;

template<typename InputType, typename ResultType>
class ServerBuffers;

template<typename InputType, typename CommonType, typename ResultType>
class RequestMembers;

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template<typename InputType, typename CommonType, typename ResultType>
class RequestManager {
	using request_t = Request<InputType, CommonType>;
	using serverBuffers_t = ServerBuffers<InputType, ResultType>;
	using requestMembers_t = RequestMembers<InputType, CommonType, ResultType>;
	using workItem_t = WorkItem<InputType, CommonType, ResultType>;
public:

	void registerRequest(std::shared_ptr<request_t> const& request);
	void removeRequest(std::shared_ptr<request_t> const& request);

	size_t queueSize() const;
	size_t containerSize() const;

	const request_t* nextRequest() noexcept;
	serverBuffers_t* buffers(std::shared_ptr<request_t> const& request);
	std::vector<workItem_t>* workItems(std::shared_ptr<request_t> const& request);
	CommonType* common(std::shared_ptr<request_t> const& request);

	bool isValid(std::shared_ptr<request_t> const& request) const;
	bool isProcessed(std::shared_ptr<request_t> const& request);

private:
	std::map<std::shared_ptr<request_t> , requestMembers_t> _requestContainer;

};

}  // namespace as

#endif /* SRC_REQUESTMANAGER_H_ */
