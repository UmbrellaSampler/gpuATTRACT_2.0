/*
 * DeviceNLGrid.h
 *
 *  Created on: Jul 21, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICENLGRID_H_
#define SRC_DEVICENLGRID_H_

#ifdef CUDA

#include "cuda_runtime.h"
#include "nativeTypes.h"
#include "Grid.h"

namespace as {

template<typename REAL>
class DeviceNLGrid {
	// Check if REAL is of floating-point type
	using real_t = typename Grid<REAL>::real_t;
	using real3_t = typename Grid<REAL>::real3_t;
	using size3_t = typename Grid<REAL>::size3_t;


public:
	struct Desc {
		size3_t dimN;				/** number of grid elements in each dimension */
		real_t dVox;				/** voxel distance */
		real_t dVox_inv;			/** inverse voxel distance */
		real_t dPlateau2;     		/** Plateau distance */
		real3_t minDim;				/** lower bound of grid dimensions */
		real3_t maxDim;				/** upper bound of grid dimensions */
		cudaTextureObject_t tex; 	/** texture object */
		unsigned*
			neighborList = nullptr; /** device NL pointer*/

		cudaTextureObject_t*
			texArrayLin = nullptr; 	/** texture objects for built-in interpolation  */
		cudaTextureObject_t*
			texArrayPt = nullptr; 	/** texture objects for manual interpolation  */
	};

	struct HostResc {
		/** Host pointers / objects*/
		cudaTextureObject_t tex;
		cudaArray* cuArray; /** cudaArray */
		unsigned* d_neighborList = nullptr;
	};

	Desc desc;
	HostResc hostResc;

};
}

#endif //CUDA



#endif /* SRC_DEVICENLGRID_H_ */
