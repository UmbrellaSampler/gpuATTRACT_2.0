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

namespace as {

template<typename BufferType>
BufferType* HostPinnedAllocator<BufferType>::allocate(size_t size) {
	BufferType* buffer;
	CUDA_CHECK(cudaMallocHost((void**)&buffer, size*sizeof(BufferType)));
	return buffer;
}

template<typename BufferType>
void HostPinnedAllocator<BufferType>::deallocate(BufferType* buffer, size_t size) {
	(void)size;
	CUDA_CHECK(cudaFreeHost(buffer));
}

} // namespace as
#endif

#endif /* SRC_HOSTPINNEDALLOCATOR_TPP_ */


