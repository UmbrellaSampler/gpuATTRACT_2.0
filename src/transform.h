/*
 * CoordTrafo.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_TRANSFORM_H_
#define SRC_TRANSFORM_H_

#include "Vec3.h"
#include "RotMat.h"
#include "MatrixFunctions.h"

#ifdef CUDA
#include "Types_6D.h"
#endif

namespace as {

template<typename REAL>
void rotate_translate(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		Vec3<REAL> const& displacement,
		Vec3<REAL> const& ang,
		unsigned const& ligSize,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)
{
	const RotMat<REAL> rotMat = euler2rotmat(ang.x, ang.y, ang.z);
	for (unsigned i = 0; i < ligSize; ++i) {
		Vec3<REAL> posAtom(x[i], y[i], z[i]);

		posAtom = rotMat*posAtom;
		posAtom += displacement;

		xTr[i] = posAtom.x;
		yTr[i] = posAtom.y;
		zTr[i] = posAtom.z;
	}
}

#ifdef CUDA

template<typename REAL>
void d_DOF2Pos(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		REAL const* x,
		REAL const* y,
		REAL const* z,
		typename Types_6D<REAL>::DOF* dofs,
		unsigned ligSize,
		unsigned numDOFs,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr);

#endif

} // namespace as




#endif /* SRC_TRANSFORM_H_ */
