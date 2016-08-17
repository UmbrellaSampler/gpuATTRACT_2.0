/*
 * trilinIntrpl.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_TRILININTRPL_H_
#define SRC_TRILININTRPL_H_

#include <type_traits>
#include "VoxelOctet.h"
#include "nativeTypesWrapper.h"
#include "nativeTypesFunctions.h"
#include "nativeTypesMath.h"


//#ifndef INLINE
//#ifdef __CUDACC__
//#define INLINE __forceinline__ __device__
//#else
//#define INLINE inline
//#endif
//#endif

namespace as {

template<typename REAL>
typename TypeWrapper<REAL>::real4_t trilinearInterpolation(
		typename TypeWrapper<REAL>::real3_t const& pos,
		const VoxelOctet<REAL> &voxelOct,
		const REAL &voxelVol_inv)
{
	using real_t = typename TypeWrapper<REAL>::real_t;
	using real3_t = typename TypeWrapper<REAL>::real3_t;
	using real4_t = typename TypeWrapper<REAL>::real4_t;

	real3_t pos_m_posMin = pos - voxelOct.min;
	real3_t posMax_m_pos = voxelOct.max - pos;

	REAL tmpMax_xy = (posMax_m_pos.x) * (posMax_m_pos.y);
	REAL tmpMin_xy = (pos_m_posMin.x) * (pos_m_posMin.y);

	float4 V = voxelOct.data[0][0][0] * (tmpMax_xy * (posMax_m_pos.z))
			+ voxelOct.data[1][0][0]
					* ((pos_m_posMin.x) * (posMax_m_pos.y) * (posMax_m_pos.z))
			+ voxelOct.data[0][1][0]
					* ((posMax_m_pos.x) * (pos_m_posMin.y) * (posMax_m_pos.z))
			+ voxelOct.data[0][0][1] * (tmpMax_xy * (pos_m_posMin.z))
			+ voxelOct.data[1][0][1]
					* ((pos_m_posMin.x) * (posMax_m_pos.y) * (pos_m_posMin.z))
			+ voxelOct.data[0][1][1]
					* ((posMax_m_pos.x) * (pos_m_posMin.y) * (pos_m_posMin.z))
			+ voxelOct.data[1][1][0] * (tmpMin_xy * (posMax_m_pos.z))
			+ voxelOct.data[1][1][1] * (tmpMin_xy * (pos_m_posMin.z));

	real4_t V_d = make_real4(
				static_cast<real_t>(V.x),
				static_cast<real_t>(V.y),
				static_cast<real_t>(V.z),
				static_cast<real_t>(V.w));

	V_d = V_d * voxelVol_inv;
	return V_d;
}

}



#endif /* SRC_TRILININTRPL_H_ */
