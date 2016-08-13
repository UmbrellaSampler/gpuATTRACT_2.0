/*
 * CPUBufferAllocator.cpp
 *
 *  Created on: Mar 18, 2016
 *      Author: uwe
 */

#ifndef SRC_HOSTALLOCATOR_TPP_
#define SRC_HOSTALLOCATOR_TPP_

#include "Allocator.h"

namespace as {

template<typename BufferType>
BufferType* HostAllocator<BufferType>::allocate(size_t size) {
	return new BufferType[size];
}


template<typename BufferType>
void HostAllocator<BufferType>::deallocate(BufferType* buffer) {
	delete[] buffer;
}

} // namespace as

#endif /* SRC_HOSTALLOCATOR_TPP_ */


