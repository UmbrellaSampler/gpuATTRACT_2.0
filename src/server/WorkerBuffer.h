/*
 * WorkerBuffer.h
 *
 *  Created on: Aug 12, 2016
 *      Author: uwe
 */

#ifndef SRC_WORKERBUFFER_H_
#define SRC_WORKERBUFFER_H_

#include <vector>
#include <cassert>

namespace as {

template<typename REAL, typename ALLOC = std::allocator<REAL>>
class WorkerBuffer {
public:

	WorkerBuffer() : _buffers(0), _bufferSize(0), _bufferDim(0) {}

	WorkerBuffer(std::size_t numBuffers, std::size_t size) : _buffers(numBuffers),  _bufferSize(size), _bufferDim(0) {
		auto alloc = ALLOC();
		for (auto& buffer : _buffers) {
			buffer = alloc.allocate(size);
		}
	}

	WorkerBuffer(std::size_t numBuffers, std::size_t dimBuffer, std::size_t size) : _buffers(numBuffers),  _bufferSize(size), _bufferDim(dimBuffer) {
			auto alloc = ALLOC();
			for (auto& buffer : _buffers) {
				buffer = alloc.allocate(size);
			}
		}

	WorkerBuffer(WorkerBuffer const& obj) = delete;
	WorkerBuffer& operator=(WorkerBuffer const& rhs) = delete;

	WorkerBuffer(WorkerBuffer&& obj) {
		*this = std::move(obj);
	}
	WorkerBuffer& operator= (WorkerBuffer&& obj) {
		_buffers = std::move(obj._buffers);
		_bufferSize = std::move(obj._bufferSize);
		return *this;
	}

	~WorkerBuffer() {
		auto alloc = ALLOC();
		for (auto& buffer : _buffers) {
			alloc.deallocate(buffer,_bufferSize);
		}
	}

	REAL* get(unsigned const& i) const noexcept {
		assert(i < _buffers.size());
		return _buffers[i];
	}

	REAL* getX() const noexcept {
		return get(0);
	}

	REAL* getY() const noexcept {
		return get(1);
	}

	REAL* getZ() const noexcept {
		return get(2);
	}

	REAL* getW() const noexcept {
		return get(3);
	}

	REAL* getV() const noexcept {
		return get(4);
	}

	REAL* getX(int i) const noexcept {
		return get(0) + _bufferDim * i;
	}

	REAL* getY(int i) const noexcept {
		return get(1) + _bufferDim * i;
	}

	REAL* getZ(int i) const noexcept {
		return get(2) + _bufferDim * i;
	}

	std::size_t numBuffers() const noexcept {
		return _buffers.size();
	}

	std::size_t bufferSize() const noexcept {
		return _bufferSize;
	}

private:
	std::vector<REAL*> _buffers;
	std::size_t _bufferSize;
	std::size_t _bufferDim;
};

}  // namespace as


#endif /* SRC_WORKERBUFFER_H_ */
