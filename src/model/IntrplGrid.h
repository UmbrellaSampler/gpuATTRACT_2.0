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
#include <algorithm>

#include "Grid.h"
#include "VoxelOctet.h"
#include "macros.h"

namespace as {

template<typename REAL>
class IntrplGrid : public Grid<REAL> {
	using grid_t = Grid<REAL>;
	using typename grid_t::real_t;
	using typename grid_t::real3_t;
	using grid_t::_pos;
	using grid_t::_dVox;
	using grid_t::_dimN;
	using typename grid_t::size3_t;
public:

	struct Desc;

	IntrplGrid(Desc desc);
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

	virtual void setPos(real3_t pos) noexcept {
		grid_t::setPos(pos);
		_maxDim.x = _pos.x + (_dimN.x - 1) * _dVox;
		_maxDim.y = _pos.y + (_dimN.y - 1) * _dVox;
		_maxDim.z = _pos.z + (_dimN.z - 1) * _dVox;
	}

	virtual void translate(real3_t displ) noexcept {
		auto pos = _pos + displ;
		setPos(pos);
	}

	/*
	 ** @brief: return the min grid indices according to the position
	 */
	int3 getIndex(real3_t const& idx) const noexcept {
		/* in cases where x is placed exactly at the boundary floor does not evaluate to "dimSize" - 2
		 * --> MIN (...)*/
		int3 idxOut;
		idxOut.x = std::min(static_cast<unsigned>(floor((idx.x - _pos.x)*_dVox_inv)), _dimN.x - 2);
		idxOut.y = std::min(static_cast<unsigned>(floor((idx.y - _pos.y)*_dVox_inv)), _dimN.y - 2);
		idxOut.z = std::min(static_cast<unsigned>(floor((idx.z - _pos.z)*_dVox_inv)), _dimN.z - 2);
		return idxOut;
	}

	/*
	 ** @brief: check out of bounds indexing
	 */
	bool outOfBounds_byIndex(int3 const& idx) const noexcept {
		return ((idx.x < 0 || (idx.x >= static_cast<int>(_dimN.x)-1)) || (idx.y < 0 || (idx.y >= static_cast<int>(_dimN.y)-1)) || (idx.z < 0 || (idx.z >= static_cast<int>(_dimN.z)-1)));
	}
	
	bool outOfBounds_byPos(real3_t const& pos) const noexcept {
		return (( (pos.x < minDim().x) || (pos.x > maxDim().x) ) ||
				( (pos.y < minDim().y) || (pos.y > maxDim().y) ) ||
				( (pos.z < minDim().z) || (pos.z > maxDim().z) ) );
	}

	bool notOutOfBounds_byPos(real3_t const& pos) const noexcept {
		return (( (pos.x >= minDim().x) && (pos.x <= maxDim().x) ) &&
				( (pos.y >= minDim().y) && (pos.y <= maxDim().y) ) &&
				( (pos.z >= minDim().z) && (pos.z <= maxDim().z) ) );
	}

	/*
	 ** @brief: returns the voxelOct for host interpolation given a set of coordinates
	 ** This method depends highly on the implementation of the Grid. Therefore this
	 ** method is part of the Grid.
	 */
	VoxelOctet<real_t> getVoxelByIndex(int3 const& idx, uint const& type) const noexcept {

		assert(idx.x >= 0 && idx.x < static_cast<int>(_dimN.x)-1);
		assert(idx.y >= 0 && idx.y < static_cast<int>(_dimN.y)-1);
		assert(idx.z >= 0 && idx.z < static_cast<int>(_dimN.z)-1);

		VoxelOctet<real_t> voxelOct;
		// compute absolute position of vertices
		voxelOct.min.x = idx.x * _dVox + _pos.x;
		voxelOct.min.y = idx.y * _dVox + _pos.y;
		voxelOct.min.z = idx.z * _dVox + _pos.z;
		voxelOct.max.x = voxelOct.min.x + _dVox;
		voxelOct.max.y = voxelOct.min.y + _dVox;
		voxelOct.max.z = voxelOct.min.z + _dVox;

		// fetch data from the grid
		assert(type < 99);
		uint idxGrid = type*_dimN.z*_dimN.y*_dimN.x;

		idxGrid += _dimN.x*(idx.z*_dimN.y + idx.y) + idx.x;
		voxelOct.data[0][0][0] = _grid[idxGrid];
		voxelOct.data[1][0][0] = _grid[idxGrid + 1];
		voxelOct.data[0][1][0] = _grid[idxGrid + _dimN.x];
		voxelOct.data[1][1][0] = _grid[idxGrid + _dimN.x + 1];

		idxGrid += _dimN.x*_dimN.y;
		voxelOct.data[0][0][1] = _grid[idxGrid];
		voxelOct.data[1][0][1] = _grid[idxGrid + 1];
		voxelOct.data[0][1][1] = _grid[idxGrid + _dimN.x];
		voxelOct.data[1][1][1] = _grid[idxGrid + _dimN.x + 1];

		return voxelOct;

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

	struct Desc {

		bool typemask[99]; /** mask of atom types which are supported in that grid */
		unsigned numGrids;	/** number of grids. The last Grid is the charge grid */
		size3_t dimN;		/** number of grid points in each dimensions */

		real_t gridSpacing; 	/** grid spacing */
		real3_t posMin;		/** lower coordinate bounds == position of grid[0] */


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
