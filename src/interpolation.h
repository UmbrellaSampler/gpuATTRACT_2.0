/*
 * interpolation.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_INTERPOLATION_H_
#define SRC_INTERPOLATION_H_

#include "Protein.h"
#include "IntrplGrid.h"
#include "nativeTypesFunctions.h"
#include "nativeTypesMath.h"
#include "VoxelOctet.h"
#include "trilinIntrpl.h"
#include <cmath>

namespace as {


template<typename REAL>
typename TypeWrapper<REAL>::real4_t interpolate(
		IntrplGrid<REAL> const* grid, typename TypeWrapper<REAL>::real3_t const& pos,
		unsigned const& type, REAL const& charge)
{
	using real4_t = typename TypeWrapper<REAL>::real4_t;

	int3 idx = grid->getIndex(pos);
	VoxelOctet<REAL> voxel = grid->getVoxelByIndex(idx, type);
	real4_t pot = trilinearInterpolation(pos, voxel,grid->voxelVol_inv());
	/* El.stat. - Forces/Energy */
	if (std::fabs(charge) > static_cast<REAL>(0.001)) {
		voxel = grid->getVoxelByIndex(idx, 0);
		real4_t pot_el = trilinearInterpolation(pos, voxel,
				grid->voxelVol_inv());
		pot_el = pot_el * charge;
		pot = pot + pot_el;
	}

	return pot;
}


template<typename REAL>
void potForce(
		IntrplGrid<REAL> const* innerGrid,
		IntrplGrid<REAL> const* outerGrid,
		const Protein<REAL>* prot,
		const REAL* LigPosX,
		const REAL* LigPosY,
		const REAL* LigPosZ,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{
	using real3_t = typename TypeWrapper<REAL>::real3_t;
	using real4_t = typename TypeWrapper<REAL>::real4_t;

	const unsigned numAtoms = prot->numAtoms();
	/* loop over all elements in LigPos/output */
	for (unsigned i = 0; i < numAtoms; ++i) {
		unsigned const& type = prot->mappedType()[i];
		if (type == 0)
			continue;
		real3_t pos = make_real3(LigPosX[i], LigPosY[i], LigPosZ[i]);
		const REAL& charge = prot->charge()[i];

		const REAL zero = static_cast<REAL>(0);
		real4_t pot = make_real4(zero, zero, zero, zero); //interpolated value for Van-der-Waals & electorstatic interactions
		if (innerGrid->outOfBounds_byPos(pos)) {
			if (!outerGrid->outOfBounds_byPos(pos)) {
				pot = interpolate(outerGrid, pos, type, charge);
			}
		} else {
			pot = interpolate(innerGrid, pos, type, charge);
		}

		data_out_x[i] = pot.x;
		data_out_y[i] = pot.y;
		data_out_z[i] = pot.z;
		data_out_E[i] = pot.w;
	}

	return;
}

}



#endif /* SRC_INTERPOLATION_H_ */
