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
