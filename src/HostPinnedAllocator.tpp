/*
 * HostPinnedBufferAllocator.cpp
 *
 *  Created on: Mar 18, 2016
 *      Author: uwe
 */

#ifndef SRC_HOSTPINNEDALLOCATOR_TPP_
#define SRC_HOSTPINNEDALLOCATOR_TPP_

#ifdef CUDA

#include "Allocator.h"
#include "cuda_runtime.h"
#include "macros.h"

using namespace as;

template<typename BufferType>
BufferType* HostPinnedAllocator<BufferType>::allocateBuffer(size_t size) {
	BufferType* buffer;
	CUDA_CHECK(cudaMallocHost((void**)&buffer, size*sizeof(BufferType)));
	return buffer;
}

template<typename BufferType>
void HostPinnedAllocator<BufferType>::freeBuffer(BufferType* buffer) {
	CUDA_CHECK(cudaFreeHost(buffer));
}

#endif

#endif /* SRC_HOSTPINNEDALLOCATOR_TPP_ */


