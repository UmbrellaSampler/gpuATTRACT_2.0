/*
 * CPUBufferAllocator.cpp
 *
 *  Created on: Mar 18, 2016
 *      Author: uwe
 */

#ifndef SRC_HOSTALLOCATOR_TPP_
#define SRC_HOSTALLOCATOR_TPP_

#include "Allocator.h"

using namespace as;

template<typename BufferType>
BufferType* HostAllocator<BufferType>::allocateBuffer(size_t size) {
	return new BufferType[size];
}


template<typename BufferType>
void HostAllocator<BufferType>::freeBuffer(BufferType* buffer) {
	delete[] buffer;
}

#endif /* SRC_HOSTALLOCATOR_TPP_ */


