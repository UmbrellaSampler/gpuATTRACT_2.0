/*
 * Service.h
 *
 *  Created on: Dec 30, 2015
 *      Author: uwe
 */

#ifndef SERVICE_H_
#define SERVICE_H_

#include <memory>

namespace as {

template<typename InputType, typename CommonType, typename ResultType>
class WorkItem;

template<typename BufferType>
class Allocator;

template<typename InputType, typename CommonType, typename ResultType>
class Service {
public:
	virtual ~Service() {}

	using input_t = InputType;
	using common_t = CommonType;
	using result_t = ResultType;

	using workItem_t = WorkItem<input_t, common_t, result_t>;

	void setInputAllocator(std::shared_ptr<Allocator<input_t>> allocator) {
		_inputAllocator = allocator;
	}

	void setResultAllocator(std::shared_ptr<Allocator<result_t>> allocator) {
		_resultAllocator = allocator;
	}

	Allocator<input_t>* getInputAllocator() {
		return _inputAllocator.get();
	}

	Allocator<result_t>* getResultAllocator() {
		return _resultAllocator.get();
	}

	virtual std::function<bool(workItem_t*)> createItemProcessor() = 0;
	virtual void initAllocators() = 0;

private:

	std::shared_ptr<Allocator<input_t>> _inputAllocator;
	std::shared_ptr<Allocator<result_t>> _resultAllocator;
};

}  // namespace as

#endif /* SERVICE_H_ */
