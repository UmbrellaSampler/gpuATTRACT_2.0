/*
 * ServerRequest.h
 *
 *  Created on: Mar 21, 2016
 *      Author: uwe
 */

#ifndef SRC_REQUESTMEMBERS_H_
#define SRC_REQUESTMEMBERS_H_

#include <vector>

namespace as {

template<typename InputType, typename ResultType>
class ServerBuffers;

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template<typename InputType, typename CommonType, typename ResultType>
class RequestMembers {
	using workItem_t = WorkItem<InputType, CommonType, ResultType>;
	using serverBuffers_t = ServerBuffers<InputType, ResultType>;
public:

	std::vector<workItem_t>* items() {
		return &_items;
	}

	serverBuffers_t* serverBuffers() {
		return &_serverBuffers;
	}

	CommonType* common() {
		return &_common;
	}

private:
	std::vector<workItem_t> _items;
	serverBuffers_t _serverBuffers;
	CommonType _common;
};

}  // namespace as



#endif /* SRC_REQUESTMEMBERS_H_ */
