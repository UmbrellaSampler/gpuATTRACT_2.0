/*
 * VoxelOctet.h
 *
 *  Created on: May 22, 2016
 *      Author: uwe
 */

#ifndef SRC_VOXELOCTET_H_
#define SRC_VOXELOCTET_H_

#include "nativeTypesWrapper.h"

namespace as {

/*
 ** @brief: Represents the volume spanned by 8 voxels.
 ** This is the basic entity interpolation is performed on.
 */

template<typename REAL>
struct VoxelOctet {
	using real3_t = typename TypeWrapper<REAL>::real3_t;

	float4 data[2][2][2];	/** function values at the voxels */
	real3_t min;				/** lower bound of the voxel coordinates */
	real3_t max;				/** upper bound of the voxel coordinates */
};

}  // namespace as

#endif /* SRC_VOXELOCTET_H_ */
