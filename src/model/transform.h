/*
 * CoordTrafo.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_TRANSFORM_H_
#define SRC_TRANSFORM_H_

#include "matrixFunctions.h"
#include "Vec3.h"
#include "RotMat.h"

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


template<typename REAL>
void mode_translate(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		unsigned const& atomsize,
		unsigned const& numModes,
		REAL const* dlig,
		REAL const* xModes,
		REAL const* yModes,
		REAL const* zModes,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)
{
	for (unsigned i = 0; i < atomSize; ++i) {
		Vec3<REAL> posAtom(x[i], y[i], z[i]);
		for(int mode=0;mode<numModes;mode++){
			posAtom.x+=dlig[mode]*xModes[i*numModes+mode];
			posAtom.y+=dlig[mode]*yModes[i*numModes+mode];
			posAtom.z+=dlig[mode]*zModes[i*numModes+mode];
		}
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
		DOF_6D<REAL>* dofs,
		unsigned ligSize,
		unsigned numDOFs,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr);

#endif

} // namespace as




#endif /* SRC_TRANSFORM_H_ */
