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

#ifndef INTRPLGRID_H_
#define INTRPLGRID_H_

#include <cassert>
#include <cmath>

#include "Grid.h"
#include "VoxelOctet.h"
#include "macros.h"

namespace as {

template<typename REAL>
class IntrplGrid;

/* Forward declaration of nested class template */
template<typename REAL>
template<typename _REAL>
struct IntrplGrid<REAL>::Desc<_REAL>;

template<typename REAL>
class IntrplGrid : public Grid<REAL> {
	using grid_t = Grid<real_t>;
public:
	using desc_t = Desc<real_t>;

	IntrplGrid(desc_t desc);
	virtual ~IntrplGrid();

	real_t voxelVol() const noexcept {
		return _voxelVol;
	}

	real_t voxelVol_inv() const noexcept {
		return _voxelVol_inv;
	}

	unsigned numTypes() const noexcept {
		return _numGrids;
	}

	real_t dVox_inv() const noexcept {
		return _dVox_inv;
	}

	real3_t minDim() const noexcept {
		return _pos;
	}

	real3_t maxDim() const noexcept {
		return _maxDim;
	}

	/*
	 ** @brief: returns a pointer to the respective grid.
	 */
	float4* getHostGridPtr(unsigned i) const noexcept {
		return _grid + i*_dimN.x*_dimN.y*_dimN.z;
	}

	void setPos(real3_t pos) noexcept {
		grid_t::setPos(pos);
		_maxDim[0] = _pos.x + (_dimN.x - 1) * _dVox;
		_maxDim[1] = _pos.y + (_dimN.y - 1) * _dVox;
		_maxDim[2] = _pos.z + (_dimN.z - 1) * _dVox;
	}

	/*
	 ** @brief: return the min grid indices according to the position
	 */
	void getIndex(const real_t &x, const real_t &y, const real_t &z, int &idxX, int &idxY, int &idxZ) const noexcept {
		/* in cases where x is place exactly at the boundary floor does not evaluate to "dimSize" - 2
		 * --> MIN (...)*/
		idxX = std::min(floor((x - _pos.x)*_dVox_inv), _dimN.x - 2);
		idxY = std::min(floor((y - _pos.y)*_dVox_inv), _dimN.y - 2);
		idxZ = std::min(floor((z - _pos.z)*_dVox_inv), _dimN.z - 2);
	}

	/*
	 ** @brief: check out of bounds indexing
	 */
	bool outOfBounds_byIndex(const int &idxX, const int &idxY, const int &idxZ) const noexcept {
		return ((idxX < 0 || (idxX >= static_cast<int>(_dimN.x)-1)) || (idxY < 0 || (idxY >= static_cast<int>(_dimN.y)-1)) || (idxZ < 0 || (idxZ >= static_cast<int>(_dimN.z)-1)));
	}
	
	bool outOfBounds_byPos(const real_t &x, const real_t &y, const real_t &z) const noexcept {
		return (( (x < minDim().x) || (x > maxDim().x) ) ||
				( (y < minDim().y) || (y > maxDim().y) ) ||
				( (z < minDim().z) || (z > maxDim().z) ) );
	}

	bool notOutOfBounds_byPos(const real_t &x, const real_t &y, const real_t &z) const noexcept {
		return (( (x >= minDim().x) && (x <= maxDim().x) ) &&
				( (y >= minDim().y) && (y <= maxDim().y) ) &&
				( (z >= minDim().z) && (z <= maxDim().z) ) );
	}

	/*
	 ** @brief: returns the voxelOct for host interpolation given a set of coordinates
	 ** This method depends highly on the implementation of the Grid. Therefore this
	 ** method is part of the Grid.
	 */
	void host_getVoxelByIndex(const int &idxX, const int &idxY, const int &idxZ, const uint &type, VoxelOctet<real_t> &voxelOct) const noexcept {

		assert(idxX >= 0 && idxX < static_cast<int>(_dimN.x)-1);
		assert(idxY >= 0 && idxY < static_cast<int>(_dimN.y)-1);
		assert(idxZ >= 0 && idxZ < static_cast<int>(_dimN.z)-1);

		// compute absolute position of vertices
		voxelOct.min.x = idxX * _dVox + _pos.x;
		voxelOct.min.y = idxY * _dVox + _pos.y;
		voxelOct.min.z = idxZ * _dVox + _pos.z;
		voxelOct.max.x = voxelOct.min.x + _dVox;
		voxelOct.max.y = voxelOct.min.y + _dVox;
		voxelOct.max.z = voxelOct.min.z + _dVox;

		// fetch data from the grid
		assert(type < 99);
		uint idx = type*_dimN.z*_dimN.y*_dimN.x;

		idx += _dimN.x*(idxZ*_dimN.y + idxY) + idxX;
		voxelOct.data[0][0][0] = _grid[idx];
		voxelOct.data[1][0][0] = _grid[idx + 1];
		voxelOct.data[0][1][0] = _grid[idx + _dimN.x];
		voxelOct.data[1][1][0] = _grid[idx + _dimN.x + 1];

		idx += _dimN.x*_dimN.y;

		voxelOct.data[0][0][1] = _grid[idx];
		voxelOct.data[1][0][1] = _grid[idx + 1];
		voxelOct.data[0][1][1] = _grid[idx + _dimN.x];
		voxelOct.data[1][1][1] = _grid[idx + _dimN.x + 1];

	}

	void freeHost() noexcept;

private:

	float4* _grid; 			/** data of host grids */

	unsigned _numGrids; 	/** number of grids */

	real3_t _maxDim;		/** max. grid dimension */
	real_t _voxelVol;			/** voxel volume */
	real_t _voxelVol_inv; 	/** pre-computed inverse of the volume */
	real_t _dVox_inv; 		/** pre-computed inverse of the voxel distance */

public:
	/*
	 ** @brief: Descriptions of a gradient-energy-grid
	 */

	template<typename _REAL>
	struct Desc {
		using real_t = typename Grid<_REAL>::real_t; // performs floating point type check

		bool typemask[99]; /** mask of atom types which are supported in that grid */
		unsigned numGrids;	/** number of grids. The last Grid is the charge grid */
		unsigned width;		/** number of elements along x */
		unsigned height;	/** number of elements along y */
		unsigned depth;		/** number of elements along z */

		real_t gridSpacing; 	/** grid spacing */
		real_t posMin[3];		/** lower coordinate bounds == position of grid[0] */


		/*
		 ** @data layout:
		 ** grid
		 ** size: numGrids*width*height*depth
		 ** {grid_0,grid_1,...,grid_(K-1)} with each grid defined as slices of 2D grids:
		 ** {s0,...,s(Nz-1)} with each slice defined like
		 ** {e00,...,e(Nx-1)0, e01,...,e(Nx-1)1,...,e(Nx-1)(Ny-1)}
		 ** K: number of grids
		 ** Nz: depth
		 ** Ny: height
		 ** Nx: width
		 */
		float4* grid; /** grid data. The first 3 entries of float4 contain the forces.
						  The last one contains the energy */
	};

};

} // namespace

#endif /* INTRPLGRID_H_ */
