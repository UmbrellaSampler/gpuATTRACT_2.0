/*
 * DeviceIntrplGrid.h
 *
 *  Created on: Jul 21, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEINTRPLGRID_H_
#define SRC_DEVICEINTRPLGRID_H_


#ifdef CUDA

#include "cuda_runtime.h"
#include "nativeTypes.h"
#include "Grid.h"

namespace as {

template<typename REAL>
class DeviceIntrplGrid {
	// Check if REAL is of floating-point type
	using real_t = typename Grid<REAL>::real_t;
	using real3_t = typename Grid<REAL>::real3_t;
	using size3_t = typename Grid<REAL>::size3_t;

public:
	struct Desc {
		size3_t dimN;				/** number of grid elements in each dimension */
		REAL dVox;				/** voxel distance */
		REAL dVox_inv;			/** inverse voxel distance */
		REAL voxelVol_inv;		/** inverse voxel volume */
		real3_t minDim;			/** lower bound of grid dimensions */
		real3_t maxDim;			/** upper bound of grid dimensions */
		cudaTextureObject_t*
			texArrayLin = nullptr; 		/** texture objects for built-in interpolation  */
		cudaTextureObject_t*
			texArrayPt = nullptr; 		/** texture objects for manual interpolation  */
	};

	struct HostResc {
		unsigned numArrays; 				/** number of cuArrays */
		/** Host pointers */
		cudaArray** cuArrayPtr = nullptr; 	/** Array of cudaArrays (host pointer)*/
		cudaTextureObject_t* h_texArrayLin = nullptr;
		cudaTextureObject_t* h_texArrayPt = nullptr;

		/** Device pointers */
		cudaTextureObject_t* d_texArrayLin = nullptr;
		cudaTextureObject_t* d_texArrayPt = nullptr;
	};

	Desc desc;
	HostResc hostResc;

};

template<typename REAL>
using d_IntrlpGrid = typename DeviceIntrplGrid<REAL>::Desc;

}

#endif //CUDA

#endif /* SRC_DEVICEINTRPLGRID_H_ */
