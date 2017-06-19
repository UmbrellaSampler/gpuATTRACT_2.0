/*
 * BufferManagement.h
 *
 *  Created on: Mar 20, 2016
 *      Author: uwe
 */

#ifndef SRC_BUFFERMANAGER_H_
#define SRC_BUFFERMANAGER_H_

#include <queue>
#include <vector>
#include <map>
#include "macros.h"

namespace as {

template<typename BufferType>
class Allocator;

template<typename BufferType>
class SingleBufferManager {
public:

	SingleBufferManager();
	~SingleBufferManager();
	void setAllocator(Allocator<BufferType>* allocator) noexcept {
		ASSERT(allocator != nullptr);
		_allocator = allocator;
	}

	/* get at least a buffer of size bufferSize */
	BufferType* getBuffer(size_t bufferSize);
	void returnBuffer(BufferType* buffer);

private:

	BufferType* allocate(size_t bufferSize);
	void deallocate(BufferType*);

	bool bufferAvailable() noexcept {
		return _bufferQueue.size() > 0;
	}

	size_t getBufferSize(BufferType* buffer) const {
		return _bufferSizeMap.at(buffer);
	}

	bool isOwned(BufferType* buffer) const {
		return _bufferSizeMap.find(buffer) != _bufferSizeMap.end();
	}

	Allocator<BufferType>* _allocator;

	std::queue<BufferType*> _bufferQueue;
	std::map<BufferType*, size_t> _bufferSizeMap;

};


template<typename InputType, typename ResultType>
class BufferManager {
public:
	void setInputAllocator(Allocator<InputType>* allocator) {
		_inputBuffMng.setAllocator(allocator);
	}

	void setResultAllocator(Allocator<ResultType>* allocator) {
		_resultBuffMng.setAllocator(allocator);
	}

	void setBuffer(size_t bufferSize) noexcept {
		_inputBuffMng.setBufferSize(bufferSize);
		_resultBuffMng.setBufferSize(bufferSize);
	}

	InputType* getInputBuffer(size_t bufferSize) noexcept {
		return _inputBuffMng.getBuffer(bufferSize);
	}

	ResultType* getResultBuffer(size_t bufferSize) noexcept {
		return _resultBuffMng.getBuffer(bufferSize);
	}

	void returnInputBuffer(InputType* buffer) noexcept {
		_inputBuffMng.returnBuffer(buffer);
	}

	void returnResultBuffer(ResultType* buffer) noexcept {
		_resultBuffMng.returnBuffer(buffer);
	}

private:
	SingleBufferManager<InputType> _inputBuffMng;
	SingleBufferManager<ResultType> _resultBuffMng;
};

}

#endif /* SRC_BUFFERMANAGER_H_ */
