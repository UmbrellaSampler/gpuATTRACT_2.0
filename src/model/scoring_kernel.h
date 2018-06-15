/*
 * scoring_kernel.h
 *
 *  Created on: Jun 12, 2018
 *      Author: glenn
 */

#ifndef SCORING_KERNEL_H_
#define SCORING_KERNEL_H_
#include "transform_modes.h"
#include "Protein.h"
#include "IntrplGrid.h"
#include "nativeTypesFunctions.h"
#include "nativeTypesMath.h"
#include "VoxelOctet.h"
#include "trilinIntrpl.h"
#include <cmath>

#ifdef CUDA
#include "nativeTypesWrapper.h"
#include "DeviceIntrplGrid.h"
#include "DeviceProtein.h"
#include "macros.h"
#endif
namespace as{

template<typename REAL, typename DOF_T>
 void d_score(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		d_Protein<REAL> const  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		double const radius_cutoff,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E);

}
#endif /* SCORING_KERNEL_H_ */
