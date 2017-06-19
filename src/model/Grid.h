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

#ifndef GRID_H_
#define GRID_H_


#include <type_traits>
#include "nativeTypesWrapper.h"
#include "nativeTypesMath.h"

namespace as{

template<typename REAL>
class Grid {
public:
	// Check if REAL is of floating-point type
	using real_t = typename TypeWrapper<REAL>::real_t;
	using real3_t = typename TypeWrapper<REAL>::real3_t;
	using size3_t = uint3;
	using grid_t = Grid<real_t>;

	Grid(size3_t dimN,
			real3_t pos, real_t dVox) :
		_dimN(dimN),
		_pos(pos),
		_dVox(dVox)
	{}

	virtual ~Grid() {}

	real_t dVox() const noexcept {
		return _dVox;
	}

	real3_t pos() const noexcept {
		return _pos;
	}

	size3_t dimN() const noexcept {
		return _dimN;
	}

	void setPos(real3_t pos) noexcept {
		_pos = pos;
	}

	void setNx(unsigned N) noexcept {
		_dimN.x = N;
	}

	void setNy(unsigned N) noexcept {
		_dimN.y = N;
	}

	void setNz(unsigned N) noexcept {
		_dimN.z = N;
	}

	void setDim(size3_t dim) noexcept {
		_dimN = dim;
	}

	virtual void translate(real3_t displ) noexcept {
		_pos = _pos + displ;
	}

protected:
	size3_t _dimN; /** Number of elements per dimension */

	real3_t _pos;  	/** lower bound of grid coordinates */

	real_t _dVox;		/** voxel distance; grid spacing */

};

} // namespace



#endif /* GRID_H_ */
