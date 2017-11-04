/*
 * transform.cu
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include "transform.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"

namespace as {

template<typename REAL>
__global__ void d_DOF2Pos(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		DOF_6D<REAL>* dofs,
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

		unsigned atomIdx = idx % numAtoms;
		Vec3<REAL> posAtom(x[atomIdx], y[atomIdx], z[atomIdx]);

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
		DOF_6D<REAL>* dofs,
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
		DOF_6D<float>* dofs,
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
		DOF_6D<double>* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		double* xTr,
		double* yTr,
		double* zTr);

}  // namespace as


