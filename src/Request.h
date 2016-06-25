/*
 * Request.h
 *
 *  Created on: Dec 23, 2015
 *      Author: uwe
 */

#ifndef REQUEST_H_
#define REQUEST_H_

#include "RequestArray.h"

namespace as {

template <typename InputType, typename CommonType>
class Request {
public:

	Request() = default;
	
	Request(InputType* data, size_t size, CommonType common):
		_array(data, size),
		_common(common)
	{}
	
	Request(Request const&) = default;
	Request& operator= (Request const&) = default;

	size_t size() const noexcept {
		return _array.size();
	}

	InputType* data() const noexcept {
		return _array.data();
	}

	CommonType common() const noexcept {
		return _common;
	}

	void setDataAndSize(InputType* data, size_t size) noexcept {
		_array.setDataAndSize(data, size);
	}

	void setCommon(CommonType const& common) noexcept {
		_common = common;
	}

	void reset() {
		_array = std::move(RequestArray<InputType>());
		_common = std::move(CommonType());
	}

private:
	RequestArray<InputType> _array;
	CommonType _common;
};

}  // namespace as


#endif /* REQUEST_H_ */
