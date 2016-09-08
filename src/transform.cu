/*
 * transform.cu
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include "transform.h"

#include "Vec3.h"
#include "RotMat.h"
#include "Types_6D.h"
#include "MatrixFunctions.h"
#include "macros.h"

namespace as {

template<typename REAL>
__global__ void d_DOF2Pos(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		typename Types_6D<REAL>::DOF* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)

{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < numAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		Vec3<REAL> const& displacement = dof.pos;
		Vec3<REAL> const& ang = dof.ang;

		Vec3<REAL> posAtom(x[idx], y[idx], z[idx]);

		const RotMat<REAL> rotMat = euler2rotmat(ang.x, ang.y, ang.z);

		posAtom = rotMat*posAtom;
		posAtom += displacement;

		xTr[idx] = posAtom.x;
		yTr[idx] = posAtom.y;
		zTr[idx] = posAtom.z;

	}
}

template<typename REAL>
void d_DOF2Pos(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		REAL const* x,
		REAL const* y,
		REAL const* z,
		typename Types_6D<REAL>::DOF* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)
{
	cudaVerifyKernel((
			d_DOF2Pos<<<gridSize, blockSize, 0, stream>>> (
				x,
				y,
				z,
				dofs,
				numAtoms,
				numDOFs,
				xTr,
				yTr,
				zTr
			)
		));
}

template
void d_DOF2Pos<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		float const* x,
		float const* y,
		float const* z,
		typename Types_6D<float>::DOF* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		float* xTr,
		float* yTr,
		float* zTr);

template
void d_DOF2Pos<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		double const* x,
		double const* y,
		double const* z,
		typename Types_6D<double>::DOF* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		double* xTr,
		double* yTr,
		double* zTr);

}  // namespace as


