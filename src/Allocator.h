/*
 * ServerBufferAllocator.h
 *
 *  Created on: Mar 18, 2016
 *      Author: uwe
 */

#ifndef SRC_ALLOCATOR_H_
#define SRC_ALLOCATOR_H_

#include <cstdlib>

namespace as {

template<typename BufferType>
class Allocator {
public:
	virtual ~Allocator() {}
	virtual BufferType* allocateBuffer(size_t) = 0;
	virtual void freeBuffer(BufferType*) = 0;
};

template<typename BufferType>
class HostAllocator : public Allocator<BufferType> {
public:
	BufferType* allocateBuffer(size_t);
	void freeBuffer(BufferType*);
};

template<typename BufferType>
class HostPinnedAllocator : public Allocator<BufferType> {
public:
	BufferType* allocateBuffer(size_t);
	void freeBuffer(BufferType*);
};


}  // namespace as

#endif /* SRC_ALLOCATOR_H_ */
