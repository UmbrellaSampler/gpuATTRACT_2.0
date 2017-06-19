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
#include "WorkerManager.h"

namespace as {

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

struct FillLevel {
	unsigned id;
	unsigned level;

	bool operator< (const FillLevel& rhs) const {
		return level < rhs.level;
	}

	/** pre-increment */
	FillLevel& operator++ () {
		++level;
		return *this;
	}

	/** post-increment */
	FillLevel operator++ (int) {
		auto tmp = *this;
		++level;
		return tmp;
	}
};

/**
 * Distributes each item in items s.t. the fill levels remain at equal height.
 * Precondition: each item must be fully distributalbe to at least one worker.
 * The getWorkerIds function object returns the a set of worker ids to which the workItems a distributed.
 * This function object is passed externally and may have different behavior depending on the application.
 */

template<typename InputType, typename CommonType, typename ResultType>
void Dispatcher<InputType, CommonType, ResultType>::dispatch(std::vector<workItem_t>* items) {

	CommonType const* common = (*items)[0].common();
	std::vector<as::workerId_t> workerIds = _getWorkerIds(common, _workerMng->poolSize());
	assert(workerIds.size() > 0 && workerIds.size());
	if (workerIds.size() == 1) {
		for (auto& item : *items) {
			_workerMng->pushItemToQueue(&item, workerIds[0]);
		}
	} else {
		auto fillLevels = getFillLevelsSorted(workerIds);
		assert(fillLevels.size() > 1);
		for (auto& item : *items) {
			distributeItem(item, fillLevels);
		}
	}
}

template<typename InputType, typename CommonType, typename ResultType>
std::vector<as::FillLevel> Dispatcher<InputType, CommonType, ResultType>::getFillLevelsSorted(std::vector<workerId_t> const& workerIds) const noexcept {
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

} // namespace as

#endif /* SRC_DISPATCHER_TPP_ */
