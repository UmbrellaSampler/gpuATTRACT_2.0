/*
 * RequestArray.h
 *
 *  Created on: Dec 31, 2015
 *      Author: uwe
 */

#ifndef REQUESTARRAY_H_
#define REQUESTARRAY_H_

namespace as {

template <typename InputType>
class RequestArray {
public:
	RequestArray():
		_data(nullptr),
		_size(0)
	{}

	RequestArray(InputType* data, size_t size):
		_data(data),
		_size(size)
	{}

	RequestArray(RequestArray const&) = default;

	size_t size() const noexcept{
		return _size;
	}

	InputType* data() const noexcept {
		return _data;
	}

	void setDataAndSize(InputType* data, size_t size) noexcept {
		_data = data;
		_size = size;
	}

private:
	InputType* _data;
	size_t _size;
};

}  // namespace as



#endif /* REQUESTARRAY_H_ */
