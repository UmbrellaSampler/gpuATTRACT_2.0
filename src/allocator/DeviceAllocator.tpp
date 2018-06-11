/*
 * DeviceAllocator.tpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEALLOCATOR_TPP_
#define SRC_DEVICEALLOCATOR_TPP_

#ifdef CUDA

#include "Allocator.h"
#include "cuda_runtime.h"
#include "macros.h"

namespace as {

template<typename BufferType>
BufferType* DeviceAllocator<BufferType>::allocate(size_t size) {
	BufferType* buffer;
	CUDA_CHECK(cudaMalloc((void**)&buffer, size*sizeof(BufferType)));
	CUDA_CHECK(cudaMemset(buffer, 0, size*sizeof(BufferType)));
	return buffer;
}

template<typename BufferType>
void DeviceAllocator<BufferType>::deallocate(BufferType* buffer, size_t size) {
	(void)size;
	CUDA_CHECK(cudaFree(buffer));
}

} // namespace as
#endif




#endif /* SRC_DEVICEALLOCATOR_TPP_ */
