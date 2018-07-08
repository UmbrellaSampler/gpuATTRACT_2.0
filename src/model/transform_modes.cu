#include "transform_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"

namespace as {





template<typename REAL, typename DOF_T >
__global__ void d_DOFPos_kernel(
		d_Protein<REAL> const  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		 REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		 REAL* buffer_trafoX, REAL* buffer_trafoY, REAL* buffer_trafoZ
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int num_atoms = protein.numAtoms;

	if (idx < num_atoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / num_atoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % num_atoms;

		Vec3<REAL> posAtom(	protein.xPos[atomIdx],
							protein.yPos[atomIdx],
							protein.zPos[atomIdx]);

		deform< REAL, DOF_T>( dof, posAtom, protein,  type_protein, idx, buffer_defoX[idx], buffer_defoY[idx], buffer_defoZ[idx] );
		translate_rotate< REAL, DOF_T>( dof, posAtom, type_protein );

		buffer_trafoX[idx] = posAtom.x;
		buffer_trafoY[idx] = posAtom.y;
		buffer_trafoZ[idx] = posAtom.z;

	}
}



template<typename REAL>
__global__ void d_rotateForces(
		REAL* xForce, REAL* yForce, REAL* zForce,
		DOF_6D_Modes<REAL>* dofs,
		unsigned const numAtoms,
		unsigned const numDofs
)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

		if (idx < numAtoms*numDofs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];


		Vec3<REAL> ForceAtom( xForce[idx], yForce[idx], zForce[idx] );
		const RotMat<REAL> rotMat = euler2rotmat( dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z );

		ForceAtom = rotMat*ForceAtom;

		xForce[idx] = ForceAtom.x;
		yForce[idx] = ForceAtom.y;
		zForce[idx] = ForceAtom.z;
	}
}

template<typename REAL>
__global__ void d_rotateForces(
		unsigned const idx_protein,
		REAL* inxForce,
		REAL* inyForce,
		REAL* inzForce,
		REAL* inE,
		REAL* outxForce,
		REAL* outyForce,
		REAL* outzForce,
		REAL* outE,
		DOF_MB_Modes<REAL>* dofs,
		unsigned const numAtoms,
		unsigned const numDofs
)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

		if (idx < numAtoms*numDofs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];


		Vec3<REAL> ForceAtom( inxForce[idx], inyForce[idx], inzForce[idx] );
		const RotMat<REAL> rotMat = euler2rotmat( dof.protein[idx_protein].ang.x, dof.protein[idx_protein].ang.y, dof.protein[idx_protein].ang.z );

		ForceAtom = rotMat*ForceAtom;
		inxForce[idx] = 0;
		inyForce[idx] = 0;
		inzForce[idx] = 0;
		outxForce[idx] += ForceAtom.x;
		outyForce[idx] += ForceAtom.y;
		outzForce[idx] += ForceAtom.z;
		outE[idx] += inE[idx];
		inE[idx] = 0;
	}
}



template<typename REAL>
void d_rotateForces(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		REAL* xForce,
		REAL* yForce,
		REAL* zForce,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numDOFs
)
{
	d_rotateForces<<<gridSize, blockSize, 0, stream>>> (
			xForce,
			yForce,
			zForce,
			dofs,
			numAtoms,
			numDOFs
			);
}

template<typename REAL>
void d_rotateForces(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned const idx_protein,
		REAL* inxForce,
		REAL* inyForce,
		REAL* inzForce,
		REAL* inE,
		REAL* outxForce,
		REAL* outyForce,
		REAL* outzForce,
		REAL* outE,
		DOF_MB_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numDOFs
)
{
	d_rotateForces<<<gridSize, blockSize, 0, stream>>> (
			idx_protein,
			inxForce,
			inyForce,
			inzForce,
			inE,
			outxForce,
			outyForce,
			outzForce,
			outE,
			dofs,
			numAtoms,
			numDOFs
			);
}




template<typename REAL, typename DOF_T>
 void d_DOFPos(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL> const&  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		 REAL* buffer_trafoX, REAL* buffer_trafoY, REAL* buffer_trafoZ
		)
{
	cudaVerifyKernel((
		d_DOFPos_kernel<<<gridSize, blockSize, 0, stream>>> (
		protein,
		dofs,
		numDOFs,
		type_protein,
		buffer_defoX,  buffer_defoY,  buffer_defoZ,
		buffer_trafoX,  buffer_trafoY,  buffer_trafoZ
		))
		);
}



template
void d_DOFPos<float, DOF_6D_Modes<float>>(
		 unsigned blockSize,
				unsigned gridSize,
				const cudaStream_t &stream,
		d_Protein<float> const&  protein,
		DOF_6D_Modes<float> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
void d_DOFPos<double, DOF_6D_Modes<double>>(
		 unsigned blockSize,
				unsigned gridSize,
				const cudaStream_t &stream,
		d_Protein<double> const&  protein,
		DOF_6D_Modes<double> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
void d_DOFPos<float, DOF_6D<float>>(
		 unsigned blockSize,
				unsigned gridSize,
				const cudaStream_t &stream,
		d_Protein<float> const&  protein,
		DOF_6D<float> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
void d_DOFPos<double, DOF_6D<double>>(
		 unsigned blockSize,
				unsigned gridSize,
				const cudaStream_t &stream,
		d_Protein<double> const&  protein,
		DOF_6D<double> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);





template
void d_rotateForces<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		float* xForce,
		float* yForce,
		float* zForce,
		DOF_6D_Modes<float>* dofs,
		unsigned numAtoms,
		unsigned numDOFs);

template
void d_rotateForces<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		double* xForce,
		double* yForce,
		double* zForce,
		DOF_6D_Modes<double>* dofs,
		unsigned numAtoms,
		unsigned numDOFs);

template
void d_rotateForces<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned const idx_protein,
		float* inxForce,
		float* inyForce,
		float* inzForce,
		float* inE,
		float* outxForce,
		float* outyForce,
		float* outzForce,
		float* outE,
		DOF_MB_Modes<float>* dofs,
		unsigned numAtoms,
		unsigned numDOFs);

template
void d_rotateForces<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned const idx_protein,
		double* inxForce,
		double* inyForce,
		double* inzForce,
		double* inE,
		double* outxForce,
		double* outyForce,
		double* outzForce,
		double* outE,
		DOF_MB_Modes<double>* dofs,
		unsigned numAtoms,
		unsigned numDOFs);

}  // namespace as
