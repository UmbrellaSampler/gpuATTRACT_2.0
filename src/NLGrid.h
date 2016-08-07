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

#ifndef NLGRID_H_
#define NLGRID_H_

#include <cmath>
#include <cassert>
#include "Grid.h"

namespace as {

struct NeighbourDesc {
	unsigned numEl;	/** number of elements in the neighbourlist.
	 	 	 	 If numEl == 0 then idx is ignored  */
	unsigned idx;	/** starting index in NLGridDesc.neighbourArray */
};

template<typename REAL>
class NLGrid : public Grid<REAL> {
	using grid_t = Grid<REAL>;
	using typename grid_t::real_t;
	using typename grid_t::real3_t;
	using grid_t::_pos;
	using grid_t::_dVox;
	using grid_t::_dimN;
	using typename grid_t::size3_t;
public:

	struct Desc;

	NLGrid(Desc desc);

	virtual ~NLGrid();

	real_t dPlateau2() const noexcept {
		return _dPlateau2;
	}

	real_t dPlateau2_inv() const noexcept {
		return _dPlateau2_inv;
	}

	unsigned* neighborList() const noexcept {
		return _neighborList;
	}

	unsigned neighborListSize() const noexcept {
		return _numElInLists;
	}

	real_t dVox_inv() const noexcept {
		return _dVox_inv;
	}

	real3_t minDim() const noexcept{
		return _pos;
	}

	real3_t maxDim() const noexcept {
		return _maxDim;
	}

	NeighbourDesc* grid() const noexcept {
		return _grid;
	}

	virtual void setPos(real3_t pos) noexcept {
		grid_t::setPos(pos);
		_maxDim.x = _pos.x + (_dimN.x - 1) * _dVox;
		_maxDim.y = _pos.y + (_dimN.y - 1) * _dVox;
		_maxDim.z = _pos.z + (_dimN.z - 1) * _dVox;
	}

	/*
	 ** @brief: return the min grid indices according to the position
	 */
	void getIndex(const real_t &x, const real_t &y, const real_t &z, int &idxX, int &idxY, int &idxZ) const noexcept {

		idxX = round((x - _pos.x)*_dVox_inv);
		idxY = round((y - _pos.y)*_dVox_inv);
		idxZ = round((z - _pos.z)*_dVox_inv);
		assert(idxX >= 0 && idxX <= static_cast<int>(_dimN.x)-1);
		assert(idxY >= 0 && idxY <= static_cast<int>(_dimN.y)-1);
		assert(idxZ >= 0 && idxZ <= static_cast<int>(_dimN.z)-1);
	}

	/*
	 ** @brief: check out of bounds according to indices
	 */
	bool outOfBounds(const int &idxX, const int &idxY, const int &idxZ) const noexcept {
		return ( idxX < 0 || (idxX > static_cast<int>(_dimN.x))) ||
				(idxY < 0 || (idxY > static_cast<int>(_dimN.y))) ||
				(idxZ < 0 || (idxZ > static_cast<int>(_dimN.z)));
	}

	/*
	 ** @brief: check out of bounds according to position
	 */
	bool outOfBounds(const real_t &x, const real_t &y, const real_t &z) const noexcept {
		return ((x < minDim().x || x > maxDim().x) ||
				(y < minDim().y || y > maxDim().y) ||
				(z < minDim().z || z > maxDim().z));
	}

	 const NeighbourDesc& getNeighbourDesc(const int &idxX, const int &idxY, const int &idxZ) const noexcept {

		return _grid[_dimN.x*(idxZ*_dimN.y + idxY) + idxX];
	}

	unsigned getNeighbor(const unsigned &idx) const noexcept {
		return _neighborList[idx];
	}

	bool outOfPlateau(const real_t &dr2) const noexcept {
		return dr2 > _dPlateau2;
	}

	real_t getRatio (real_t d2) const noexcept {
		return _ratio[int(d2*10000)];
	}

	/*
	 ** @brief: Free resources on the host.
	 */
	void freeHost() noexcept;

	NeighbourDesc* _grid;	/** host grid of neighbourlists */

	unsigned _numElInLists;	/** number of total elements in Lists */
	unsigned* _neighborList;	/** contains receptor atom indices for all neighbours.
								It is basically a concatenation of subsequent neighbour lists.
								The exact sequence is the same as for the grid */

	real3_t _maxDim;			/** max. grid dimension */
	real_t _dPlateau2;		/** plateau distance squared */
	real_t _dPlateau2_inv;	/** plateau distance squared inverse*/
	real_t _dVox_inv;		/** pre-calculated inverse of the voxel distance */

	real_t *_ratio;	/** precomputed square roots of the distance squared */

public:
	/*
	 ** @brief: Descriptions of a neighbour list grid
	 */
	struct Desc {

		size3_t dimN;		/** number of grid points in each dimensions */
		real_t gridSpacing; 	/** grid spacing */
		real_t dPlateau; 	/** Plateau distance: cutOff */
		real3_t posMin;		/** lower coordinate bounds == position of grid[0] */

		/*
		 ** @data layout:
		 ** grid
		 ** size: width*height*depth
		 ** slices of 2D grids {s0,...,s(Nz-1)} with each slice defined like
		 ** {e00,...,e(Nx-1)0, e01,...,e(Nx-1)1,...,e(Nx-1)(Ny-1)}
		 ** Nz: depth
		 ** Ny: height
		 ** Nx: width
		 */
		NeighbourDesc* grid;	/** grid of neighbour descriptions */

		unsigned numEl;				/** number of total elements in neighbourArray */

		/*
		 ** @data layout:
		 ** neighbourArray
		 ** size: numEl
		 ** {list_0,...,list_(N-1)} with each list
		 ** {idx_0,...,idx_Li-1} where Li is the length of the i-th list
		 ** The sum of all Li's equals numEl
		 ** N = width*height*depth
		 */
		unsigned* neighborArray;	/** contains receptor atom indices for all neighbours.
									It is basically a concatenation of subsequent neighbour lists.
									The exact sequence is the same as for the grid */
	};
};

} // namespace

#endif /* NLGRID_H_ */
