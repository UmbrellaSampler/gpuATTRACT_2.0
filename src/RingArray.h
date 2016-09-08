/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef RINGARRAY_H_
#define RINGARRAY_H_

#include <memory>

namespace as {

template <typename T>
class RingArray {
	public:
	/* Constructor */
	RingArray(size_t size) : _ring(size), _pointer(-1) {}

	void push(const T& obj) {
		rotate();
		_ring[calcIdx(0)] = obj;
	}

	T get(const unsigned& stage) const {
		return _ring[calcIdx(stage)];
	}


	void rotate() {
		_pointer = (_pointer + 1) % _ring.size();
	}

private:

	unsigned calcIdx (const unsigned& stage) const {
		return ((int)_ring.size() -_pointer + (int)stage) % _ring.size();
	}

	std::vector<T> _ring;
	int _pointer;

};

}  // namespace AServ


#endif /* RINGARRAY_H_ */
