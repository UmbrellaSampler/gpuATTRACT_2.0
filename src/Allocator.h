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
	virtual BufferType* allocate(size_t) = 0;
	virtual void deallocate(BufferType*) = 0;
};

template<typename BufferType>
class HostAllocator : public Allocator<BufferType> {
public:
	BufferType* allocate(size_t);
	void deallocate(BufferType*);
};

template<typename BufferType>
class HostPinnedAllocator : public Allocator<BufferType> {
public:
	BufferType* allocate(size_t);
	void deallocate(BufferType*);
};


}  // namespace as

#endif /* SRC_ALLOCATOR_H_ */
