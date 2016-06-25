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

#ifndef SRC_INTRPLGRID_TPP_
#define SRC_INTRPLGRID_TPP_

#include "IntrplGrid.h"

using namespace as;

template<typename REAL>
IntrplGrid<REAL>::IntrplGrid(desc_t desc):
		Grid<REAL>::Grid(desc.width, desc.height,	desc.depth,
				real3_t(desc.posMin[0], desc.posMin[1], desc.posMin[2]),
				desc.gridSpacing),
		_grid(desc.grid),
		_numGrids(desc.numGrids),
		_voxelVol(_dVox * _dVox * _dVox),
		_voxelVol_inv(1.f/_voxelVol),
		_dVox_inv(1.f/_dVox)

{
	_maxDim[0] = _pos.x + (_dimN.x - 1) * _dVox;;
	_maxDim[1] = _pos.y + (_dimN.y - 1) * _dVox;;
	_maxDim[2] = _pos.z + (_dimN.z - 1) * _dVox;;
}

template<typename REAL>
IntrplGrid<REAL>::~IntrplGrid() {
	freeHost();
}

template<typename REAL>
void IntrplGrid<REAL>::freeHost() noexcept {
	delete[] _grid;
}

#endif


/* associated functions */
//only used for debugging

// TODO: create & move to in GridUtilities.h/.cpp

//template<typename REAL>
//void print(IntrplGrid<REAL>* desc, uint grid0, uint gridM,
//		uint comp0, uint compM,
//		uint x0, uint xM, uint y0, uint yM, uint z0, uint zM)
//{
//	// Store io flags
//	using namespace std;
//	int precisionSetting = cout.precision( );
//	ios::fmtflags flagSettings = cout.flags();
//
//	cout.setf(ios::fixed | ios::showpos | ios::showpoint);
//
//	cout << "numGrids " << desc->numTypes() << endl;
//	cout << "width " << desc->width() << endl;
//	cout << "height " << desc->height() << endl;
//	cout << "depth " << desc->depth() << endl;
//	cout << "dVox " << desc->dVox() << endl;
//	cout << "minPos " << desc->pos().x << " " << desc->pos().y << " " << desc->pos().z << endl;
//
////	cout << "Debug"<< " " <<   grid0 << " " <<  gridM << " " <<  comp0 << " " <<  compM << " " <<  x0 << " " <<  xM << " " <<  y0 << " " <<  yM << " " <<  z0 << " " <<  zM  << endl;
//
//	cout.precision(6);
//	uint xOff = 5;
//	uint yOff = 2;
//	uint width = 13;
//	uint height = 2; // > 0 !!!
//	for (uint grid = grid0; grid <= gridM; ++grid) {
//		float4* grid_ptr = desc->getHostGridPtr(grid);
//		cout << "###GRID " << grid <<"###" << endl;
//		for (uint comp = comp0; comp <= compM; ++comp) {
//			cout << "###COMP " << comp <<"###" << endl;
//			for (uint z = z0; z <= zM; ++z) {
//				cout << "---Slice" << z << "---" << endl;
//				cout << setw(xOff) << " ";
//				for (uint x = x0; x <= xM; ++x)
//					cout << setw(width)<<  x;
//				for (uint i = 0; i < yOff; ++i)
//					cout << endl;
//
//				for (uint y = y0; y <= yM; ++y) {
//					cout << setw(xOff) <<  y;
//					for (uint x = x0; x <= xM; ++x) {
//						uint idx =  desc->width()*(z*desc->height() + y) + x;
//						float* ptr = (float*)&grid_ptr[idx];
//						cout << setw(width) << ptr[comp];
//					}
//					for (uint i = 0; i < height; ++i)
//						cout << endl;
//				}
//				cout << endl;
//			}
//		}
//	}
//
//	// Restore io flags
//	cout.precision(precisionSetting);
//	cout.flags(flagSettings);
//
//}
//template<typename REAL>
//void print(const GradEnGridDesc<REAL>& desc, uint grid0, uint gridM, uint comp0, uint compM, uint x0, uint xM, uint y0, uint yM, uint z0, uint zM)
//{
//	// Store io flags
//	using namespace std;
//	int precisionSetting = cout.precision( );
//	ios::fmtflags flagSettings = cout.flags();
//
//	cout.setf(ios::fixed | ios::showpos | ios::showpoint);
//
//	cout << "numGrids " << desc.numGrids << endl;
//	cout << "width " << desc.width << endl;
//	cout << "height " << desc.height << endl;
//	cout << "depth " << desc.depth << endl;
//	cout << "dVox " << desc.gridSpacing << endl;
//	cout << "minPos " << desc.posMin[0] << " " << desc.posMin[1] << " " << desc.posMin[2] << endl;
//
////	cout << "Debug"<< " " <<   grid0 << " " <<  gridM << " " <<  comp0 << " " <<  compM << " " <<  x0 << " " <<  xM << " " <<  y0 << " " <<  yM << " " <<  z0 << " " <<  zM  << endl;
//
//	cout.precision(6);
//	uint xOff = 5;
//	uint yOff = 2;
//	uint width = 13;
//	uint height = 2; // > 0 !!!
//	for (uint grid = grid0; grid <= gridM; ++grid) {
//		float4* grid_ptr = desc.grid + grid*desc.width*desc.height*desc.depth;
//		cout << "###GRID " << grid <<"###" << endl;
//		for (uint comp = comp0; comp <= compM; ++comp) {
//			cout << "###COMP " << comp <<"###" << endl;
//			for (uint z = z0; z <= zM; ++z) {
//				cout << "---Slice" << z << "---" << endl;
//				cout << setw(xOff) << " ";
//				for (uint x = x0; x <= xM; ++x)
//					cout << setw(width)<<  x;
//				for (uint i = 0; i < yOff; ++i)
//					cout << endl;
//
//				for (uint y = y0; y <= yM; ++y) {
//					cout << setw(xOff) <<  y;
//					for (uint x = x0; x <= xM; ++x) {
//						uint idx =  desc.width*(z*desc.height + y) + x;
//						float* ptr = (float*)&grid_ptr[idx];
//						cout << setw(width) << ptr[comp];
//					}
//					for (uint i = 0; i < height; ++i)
//						cout << endl;
//				}
//				cout << endl;
//			}
//		}
//	}
//
//	// Restore io flags
//	cout.precision(precisionSetting);
//	cout.flags(flagSettings);
//
//}



