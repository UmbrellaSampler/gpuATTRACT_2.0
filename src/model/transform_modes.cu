#include "transform_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"

namespace as {


template <typename REAL, typename DOF_T>
__device__ void deform( DOF_T const& dof, Vec3<REAL> & posAtom, d_Protein<REAL> const*  protein, unsigned const idxAtom, unsigned const idxModes, unsigned const bufIdx,
		REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	unsigned const numModes = protein->numModes;
	for(int mode=0; mode < numModes; mode++){
		posAtom.x += dof.modesRec[mode] * protein->xModes[idxAtom*numModes+mode];
		posAtom.y += dof.modesRec[mode] * protein->yModes[idxAtom*numModes+mode];
		posAtom.z += dof.modesRec[mode] * protein->zModes[idxAtom*numModes+mode];
	}
		buffer_defoX[bufIdx] = posAtom.x;
		buffer_defoZ[bufIdx] = posAtom.y;
		buffer_defoZ[bufIdx] = posAtom.z;

}

template <typename REAL, typename DOF_T>
__device__ void deform( DOF_T const& dof, Vec3<REAL> & posAtom, d_Protein<REAL> const*  protein, unsigned const idxAtom, unsigned const idxModes, unsigned const bufIdx,
		REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0 )
{

}


template <typename REAL, typename DOF_T, int PROTEIN_T>
__device__ void translate_rotate( DOF_T const& dof, Vec3<REAL> & posAtom,
	typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
	if ( PROTEIN_T == 0){
		rotMat = rotMat.getInv();
	}

	posAtom = rotMat*posAtom;
	posAtom += dof._6D.pos;
}

template <typename REAL, typename DOF_T, int PROTEIN_T>
__device__ void translate_rotate( DOF_T const& dof, Vec3<REAL> & posAtom ,
	typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0 )
{
	RotMat<REAL> rotMat = euler2rotmat(dof.ang.x, dof.ang.y, dof.ang.z);
	if ( PROTEIN_T == 0){
		rotMat = rotMat.getInv();
	}
	posAtom = rotMat*posAtom;
	posAtom += dof.pos;
}

template<typename REAL, typename DOF_T, int PROTEIN_T >
__global__ void d_DOFPos_kernel(
		d_Protein<REAL> const*  protein,
		DOF_T* dofs,
		unsigned const numDOFs,
		 REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		 REAL* buffer_trafoX, REAL* buffer_trafoY, REAL* buffer_trafoZ
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int num_atoms = protein->numAtoms;

	if (idx < num_atoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / num_atoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % num_atoms;

		Vec3<REAL> posAtom(	protein->xPos[atomIdx],
							protein->xPos[atomIdx],
							protein->xPos[atomIdx]);

		deform< REAL, DOF_T>( dof, posAtom, protein, atomIdx, 0, idx, buffer_defoX, buffer_defoY, buffer_defoZ);

		translate_rotate< REAL, DOF_T, PROTEIN_T>( dof, posAtom );

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


		Vec3<REAL> ForceAtom(xForce[idx], yForce[idx], zForce[idx]);
		const RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);

		ForceAtom = rotMat*ForceAtom;

		xForce[idx] = ForceAtom.x;
		yForce[idx] = ForceAtom.y;
		zForce[idx] = ForceAtom.z;
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




template<typename REAL, typename DOF_T, int PROTEIN_T >
 void d_DOFPos(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL> const*  protein,
		DOF_T* dofs,
		unsigned const numDOFs,
		REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		 REAL* buffer_trafoX, REAL* buffer_trafoY, REAL* buffer_trafoZ
		)
{
	cudaVerifyKernel((
		d_DOFPos_kernel< REAL, DOF_T,  PROTEIN_T ><<<gridSize, blockSize, 0, stream>>> (
		blockSize,
		gridSize,
		stream,
		protein,
		dofs,
		numDOFs,
		buffer_defoX,  buffer_defoY,  buffer_defoZ,
		buffer_trafoX,  buffer_trafoY,  buffer_trafoZ
		))
		);
}


template
__global__ void d_DOFPos_kernel<float, DOF_6D_Modes<float>, 0>(

		d_Protein<float> const*  protein,
		DOF_6D_Modes<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<double, DOF_6D_Modes<double>, 0>(

		d_Protein<double> const*  protein,
		DOF_6D_Modes<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<float, DOF_6D_Modes<float>, 1>(

		d_Protein<float> const*  protein,
		DOF_6D_Modes<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<double, DOF_6D_Modes<double>, 1>(

		d_Protein<double> const*  protein,
		DOF_6D_Modes<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<float, DOF_6D<float>, 0>(

		d_Protein<float> const*  protein,
		DOF_6D<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<double, DOF_6D<double>, 0>(

		d_Protein<double> const*  protein,
		DOF_6D<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<float, DOF_6D<float>, 1>(

		d_Protein<float> const*  protein,
		DOF_6D<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
__global__ void d_DOFPos_kernel<double, DOF_6D<double>, 1>(

		d_Protein<double> const*  protein,
		DOF_6D<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

/*
template
 void d_DOFPos<float, DOF_6D_Modes<float>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D_Modes<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
 void d_DOFPos<double, DOF_6D_Modes<double>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D_Modes<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
 void d_DOFPos<float, DOF_6D_Modes<float>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D_Modes<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
 void d_DOFPos<double, DOF_6D_Modes<double>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D_Modes<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
 void d_DOFPos<float, DOF_6D<float>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
 void d_DOFPos<double, DOF_6D<double>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);

template
 void d_DOFPos<float, DOF_6D<float>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D<float>* dofs,
		unsigned const numDOFs,
		float* buffer_defoX, float* buffer_defoY, float* buffer_defoZ,
		float* buffer_trafoX, float* buffer_trafoY, float* buffer_trafoZ
		);

template
 void d_DOFPos<double, DOF_6D<double>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D<double>* dofs,
		unsigned const numDOFs,
		double* buffer_defoX, double* buffer_defoY, double* buffer_defoZ,
		double* buffer_trafoX, double* buffer_trafoY, double* buffer_trafoZ
		);


*/




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

}  // namespace as
