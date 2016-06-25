/*
 * Dispatcher.tpp
 *
 *  Created on: Mar 24, 2016
 *      Author: uwe
 */

#ifndef SRC_DISPATCHER_TPP_
#define SRC_DISPATCHER_TPP_

#include <cassert>
#include <algorithm>
#include "macros.h"
#include "Dispatcher.h"
#include "Service.h"

using namespace as;

template<typename InputType, typename CommonType, typename ResultType>
void Dispatcher<InputType, CommonType, ResultType>::run() {
	assert(_workerMng != nullptr);
	assert(_requestQueue != nullptr);
	while(true) {
		auto items = _requestQueue->waitAndPop();
		if (items == nullptr) {
			break;
		}
		assert(items->size() > 0);

		dispatch(items);

	}
}

struct as::FillLevel {
	unsigned id;
	unsigned level;

	bool operator< (const FillLevel& rhs) const {
		return level < rhs.level;
	}

	/* pre-increment */
	FillLevel& operator++ () {
		++level;
		return *this;
	}

	/* post-increment */
	FillLevel operator++ (int) {
		auto tmp = *this;
		++level;
		return tmp;
	}
};

template<typename InputType, typename CommonType, typename ResultType>
void Dispatcher<InputType, CommonType, ResultType>::dispatch(std::vector<workItem_t>* items) {

	if (_workerMng->poolSize() == 1) {
		for (auto& item : *items) {
			_workerMng->pushItemToQueue(&item, 0);
		}
	} else {
		CommonType const* common = (*items)[0].common();
		std::set<unsigned> workerIds = getWorkerIds(common);
		auto fillLevels = getFillLevelsSorted(workerIds);
		assert(fillLevels.size() > 1);
		for (auto& item : *items) {
			distributeItem(item, fillLevels);
		}
	}
}

template<typename InputType, typename CommonType, typename ResultType>
std::set<unsigned> Dispatcher<InputType, CommonType, ResultType>::getWorkerIds(CommonType const* common) const noexcept{
	// TODO: Note: this is a temporary solution
	//
	/* compiler dummy */ ASSERT(common);
	std::set<unsigned> workerIds;
	size_t size = _workerMng->poolSize();
	for (unsigned i = 0; i < size; ++i) {
		workerIds.insert(i);
	}
	// Now get Ids with dataMng!
	return workerIds;
}

template<typename InputType, typename CommonType, typename ResultType>
std::vector<as::FillLevel> Dispatcher<InputType, CommonType, ResultType>::getFillLevelsSorted(std::set<unsigned> const& workerIds) const noexcept {
	const auto size = workerIds.size();
	std::vector<FillLevel> levels(size);
	unsigned count = 0;
	for (auto const& id : workerIds) {
		FillLevel& level = levels[count];
		level.id = id;
		level.level = _workerMng->queueSize(id);
		++count;
	}
	std::sort(levels.begin(), levels.end());
	return levels;
}

template<typename InputType, typename CommonType, typename ResultType>
void Dispatcher<InputType, CommonType, ResultType>::distributeItem(workItem_t& item, std::vector<FillLevel>& levels) {
	size_t numWorkers = levels.size();
	for (size_t i = 0; i < numWorkers-1; ++i) {
		auto& currLevel = levels[i];
		auto& nextLevel = levels[i+1];
		if (currLevel < nextLevel) {
			_workerMng->pushItemToQueue(&item, currLevel.id);
			++currLevel;
			break;
		}
		if (i == numWorkers-2) {
			_workerMng->pushItemToQueue(&item, nextLevel.id);
			++nextLevel;
		}
	}
}

#endif /* SRC_DISPATCHER_TPP_ */
