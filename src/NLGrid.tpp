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

#ifndef SRC_NLGRID_TPP_
#define SRC_NLGRID_TPP_

#include <cmath>

#include "NLGrid.h"

using namespace as;

template<typename REAL>
NLGrid<REAL>::NLGrid(desc_t desc):
		Grid<REAL>::Grid(desc.width, desc.height, desc.depth,
				real3_t(desc.posMin[0], desc.posMin[1], desc.posMin[2]),
				desc.gridSpacing),
		_grid(desc.grid),
		_numElInLists(desc.numEl),
		_neighborList(desc.neighborArray),
		_dPlateau2(desc.dPlateau*desc.dPlateau),
		_dPlateau2_inv(1/_dPlateau2),
		_dVox_inv(1/_dVox)
{
	_maxDim.x = _pos.x + (_dim.x - 1) * _dVox;
	_maxDim.y = _pos.y + (_dim.y - 1) * _dVox;
	_maxDim.z = _pos.z + (_dim.z - 1) * _dVox;

	/* init ratios (taken from the original attract code) */
	int size  = int(10000*_dPlateau2);
	_ratio = new real_t[size+1];

	for (int n = 0; n <= size; n++) {
		double d2 = ((n + 0.5) / 10000);
		_ratio[n] = sqrt(d2 / _dPlateau2);
	}
}

template<typename REAL>
NLGrid<REAL>::~NLGrid() {
	freeHost();
}

template<typename REAL>
void NLGrid<REAL>::freeHost() noexcept {
	delete[] _neighborList;
	delete[] _grid;
	delete[] _ratio;
}

#endif


