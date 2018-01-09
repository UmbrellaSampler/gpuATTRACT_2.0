/*
 * interpolation.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_INTERPOLATION_MB_MODES_H_
#define SRC_INTERPOLATION_MB_MODES_H_

#include "Protein.h"
#include "IntrplGrid.h"
#include "nativeTypesFunctions.h"
#include "nativeTypesMath.h"
#include "VoxelOctet.h"
#include "trilinIntrpl.h"
#include <cmath>
#include "interpolation.h"
#ifdef CUDA
#include "DeviceIntrplGrid.h"
#include "DeviceProtein.h"
#include "Types_6D_MB_Modes.h"
#endif

namespace as {


template<typename REAL>
void potForce(
		IntrplGrid<REAL> const* innerGrid,
		IntrplGrid<REAL> const* outerGrid,
		const Protein<REAL>* prot,
		DOF_6D_MB_Modes<REAL> dof,
		unsigned int & ligIdx,
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

	const RotMat<REAL> rotMat = euler2rotmat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);

	const unsigned numAtoms = prot->numAtoms();
	/* loop over all elements in LigPos/output */
	for (unsigned i = 0; i < numAtoms; ++i) {
		unsigned const& type = prot->mappedType()[i];
		constexpr REAL zero = static_cast<REAL>(0);
		real4_t pot = make_real4(zero, zero, zero, zero); //interpolated value for Van-der-Waals & electorstatic interactions
		if (type != 0) {
			const real3_t pos = make_real3(LigPosX[i], LigPosY[i], LigPosZ[i]);
			const REAL& charge = prot->charge()[i];

			if (innerGrid->outOfBounds_byPos(pos)) {
				if (!outerGrid->outOfBounds_byPos(pos)) {
					pot = interpolate(outerGrid, pos, type, charge);
				}
			} else {
//				static int count = 0;
//				if (++count < 50) {
//					printf("%u\n", i);
//				}
				pot = interpolate(innerGrid, pos, type, charge);
			}
		}

		Vec3<REAL> force(pot.x, pot.y, pot.z);
		force = rotMat * force;

		data_out_x[i] += force.x;
		data_out_y[i] += force.y;
		data_out_z[i] += force.z;
		data_out_E[i] += pot.w;
	}

	return;
}



#ifdef CUDA

template<typename REAL>
void d_potForce (
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned int ligIdx,
		const d_IntrlpGrid<REAL>& inner,
		const d_IntrlpGrid<REAL>& outer,
		const d_Protein<REAL>& prot,
		const unsigned& numDOFs,
		const REAL* data_in_x,
		const REAL* data_in_y,
		const REAL* data_in_z,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E);

#endif

} // namespace





#endif /* SRC_INTERPOLATION_MB_MODES_H_ */
