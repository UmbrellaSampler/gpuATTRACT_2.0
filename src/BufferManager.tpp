/*
 * BufferManager.cpp
 *
 *  Created on: Mar 20, 2016
 *      Author: uwe
 */

#ifndef SRC_BUFFERMANAGER_TPP_
#define SRC_BUFFERMANAGER_TPP_

#include "BufferManager.h"
#include "macros.h"
#include <cassert>
#include "Allocator.h"

using namespace as;

template<typename BufferType>
SingleBufferManager<BufferType>::SingleBufferManager() :
	_allocator(nullptr)
{}

template<typename BufferType>
SingleBufferManager<BufferType>::~SingleBufferManager() {
	/* All buffers are returned */
	ASSERT(_bufferQueue.size() == _bufferSizeMap.size() && "Have you called 'wait(...)' for all submitted requests?");

	for (auto& pair : _bufferSizeMap) {
		delete[] pair.first;
	}
}

template<typename BufferType>
inline BufferType* SingleBufferManager<BufferType>::allocateBuffer(size_t bufferSize) {
	assert(_allocator != nullptr);
	BufferType* buffer = _allocator->allocateBuffer(bufferSize);
	assert(!isOwned(buffer));
	_bufferSizeMap[buffer] = bufferSize;
	return buffer;

}

template<typename BufferType>
inline void SingleBufferManager<BufferType>::freeBuffer(BufferType* buffer) {
	assert(isOwned(buffer));
	_bufferSizeMap.erase(buffer);
	_allocator->freeBuffer(buffer);

}


template<typename BufferType>
BufferType* SingleBufferManager<BufferType>::getBuffer(size_t bufferSize) {
	assert(bufferSize > 0);
	if (!bufferAvailable()) {
		return allocateBuffer(bufferSize);
	} else {
		BufferType* buffer = _bufferQueue.front();
		_bufferQueue.pop();
		if(getBufferSize(buffer) < bufferSize) {
			freeBuffer(buffer);
			return allocateBuffer(bufferSize);
		} else {
			return buffer;
		}
	}
}

template<typename BufferType>
void SingleBufferManager<BufferType>::returnBuffer(BufferType* buffer) {
	assert(buffer != nullptr);
	assert(isOwned(buffer));
	_bufferQueue.push(buffer);
}

#endif /* SRC_BUFFERMANAGER_TPP_ */

