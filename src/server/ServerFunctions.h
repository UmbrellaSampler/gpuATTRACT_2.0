/*
 * ServerFunctions.h
 *
 *  Created on: Dec 31, 2015
 *      Author: uwe
 */

#ifndef SERVERFUNCTIONS_H_
#define SERVERFUNCTIONS_H_

#include <vector>
#include <algorithm>
#include <cassert>
#include "WorkItem.h"
#include "ServerBuffers.h"

namespace as {

inline unsigned calcNumberItems(unsigned requestSize, unsigned itemSize) {
	return (requestSize + itemSize - 1) / itemSize;
}

template <typename InputType, typename CommonType, typename ResultType>
void
initializeWorkItems
	(std::vector<as::WorkItem<InputType, CommonType, ResultType>>& items,
			as::ServerBuffers<InputType, ResultType> const* serverBuffers,
			const size_t& requestSize,
			CommonType* common_ptr,
			const size_t& itemSize
	)
{
	const unsigned numItems = items.size();
	for (unsigned i = 0; i < numItems-1; ++i) {
		auto& item = items[i];
		item.setDataAndSize(serverBuffers->inputBuffer + i*itemSize, itemSize);
		item.setResultPtr(serverBuffers->resultBuffer + i*itemSize);
		item.setCommonPtr(common_ptr);
	}
	const unsigned shift = (numItems - 1)*itemSize;
	const unsigned itemSizeLast = requestSize - shift;
	auto& itemLast = items.back();
	itemLast.setDataAndSize(serverBuffers->inputBuffer + shift, itemSizeLast);
	itemLast.setResultPtr(serverBuffers->resultBuffer + shift);
	itemLast.setCommonPtr(common_ptr);
}

template <typename InputType, typename CommonType, typename ResultType>
std::vector<as::WorkItem<InputType, CommonType, ResultType>>
createWorkItemsFn
	(as::ServerBuffers<InputType, ResultType> const* serverBuffers,
			size_t const& requestSize,
			CommonType* common,
			size_t const& itemSize)
{
	auto numItems = calcNumberItems(requestSize, itemSize);
	std::vector<as::WorkItem<InputType, CommonType, ResultType>> items(numItems);
	initializeWorkItems(items, serverBuffers, requestSize, common, itemSize);
	return items;
}

} // namespace as

#endif /* SERVERFUNCTIONS_H_ */
