/*
 * WorkItem.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef WORKITEM_H_
#define WORKITEM_H_

#include <atomic>
#include "RequestArray.h"


namespace as {

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem {
public:
	WorkItem() :
		_common_ptr(nullptr),
		_result(nullptr),
		_processed(false)
	{}

	WorkItem(InputType* array, unsigned size, CommonType* common_ptr, ResultType* result) :
		_common_ptr(common_ptr),
		_array(array, size),
		_result(result),
		_processed(false)
	{}

	WorkItem(WorkItem const&) = delete;
	WorkItem& operator=(WorkItem const&) = delete;

	unsigned size() const noexcept {
		return _array.size();
	}

	InputType* inputBuffer() const noexcept {
		return _array.data();
	}

	CommonType* common() const noexcept {
		return _common_ptr;
	}

	ResultType* resultBuffer() const noexcept {
		return _result;
	}

	void setProcessed() noexcept {
		std::atomic_thread_fence(std::memory_order_release);
		_processed.store(true,std::memory_order_relaxed);
	}

	bool isProcessed() const noexcept {
		bool ready = _processed.load(std::memory_order_relaxed);
//		std::atomic_thread_fence(std::memory_order_acquire);
		return ready;
	}

	void setDataAndSize(InputType* data, unsigned size) noexcept {
		_array.setDataAndSize(data, size);
	}

	void setCommonPtr(CommonType* common_ptr) noexcept {
		_common_ptr = common_ptr;
	}

	void setResultPtr(ResultType* result) noexcept {
		_result = result;
	}

	bool isValid() const noexcept {
		bool validArray = (_array.data() != nullptr) && (_array.size() > 0);
		bool validCommon = _common_ptr != nullptr;
		bool validResult = _result != nullptr;
		return validArray && validCommon && validResult;
	}

private:
	RequestArray<InputType> _array;
	CommonType* _common_ptr;
	ResultType* _result;

	std::atomic<bool> _processed;
};

}  // namespace as




#endif /* WORKITEM_H_ */
