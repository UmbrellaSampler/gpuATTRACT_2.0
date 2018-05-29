#include "transform_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"
#include "Allocator.h"

namespace as {

template<typename REAL>
using Buffer = WorkerBuffer<REAL, DeviceAllocator<REAL>>;


template <typename REAL, typename DOF_T>
std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type
__device__ deform( DOF_T const& dof, Vec3<REAL> & posAtom, d_Protein<REAL> const*  protein, unsigned const idxAtom, unsigned const idxModes ) {
	unsigned const numModes = protein->numModes;
	for(int mode=0; mode < numModes; mode++){
		posAtom.x += dof.modes[idxModes][mode] * xModes[idxAtom*numModes+mode];
		posAtom.y += dof.modes[idxModes][mode] * yModes[idxAtom*numModes+mode];
		posAtom.z += dof.modes[idxModes][mode] * zModes[idxAtom*numModes+mode];
	}
}

template <typename REAL, typename DOF_T>
std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type
__device__ deform( DOF_T const& dof, Vec3<REAL> & posAtom, d_Protein<REAL> const&  protein ) {

}

template <typename REAL, typename DOF_T, int PROTEIN_T>
std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type
__device__ translate_rotate( DOF_T const& dof, Vec3<REAL> & posAtom ) {
	const RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
	posAtom = rotMat*posAtom;
	posAtom += dof._6D.pos;
}

template <typename REAL, typename DOF_T, int PROTEIN_T>
std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type
__device__ translate_rotate( DOF_T const& dof, Vec3<REAL> & posAtom ) {
	const RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
	posAtom = rotMat*posAtom;
	posAtom += dof.pos;
}

template<typename REAL, typename DOF_T, int PROTEIN_T >
__global__ void d_DOFPos(
		d_Protein<REAL> const*  protein,
		DOF_T* dofs,
		unsigned const numDOFs,
		Buffer<REAL> buffer_defo,
		Buffer<REAL> buffer_trafo
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

		Vec3<REAL> posAtom(protein.xPos[atomIdx], protein.xPos[atomIdx], protein.xPos[atomIdx]);

		deform< REAL, DOF_T>( dof, posAtom, protein, idxAtom, idxModes );

		buffer_defo.getX()[bufIdx] = posAtom.x;
		buffer_defo.getY()[bufIdx] = posAtom.y;
		buffer_defo.getZ()[bufIdx] = posAtom.z;

		translate_rotate< REAL, DOF_T, PROTEIN_T>( dof, posAtom );

		buffer_trafo.getX()[bufIdx] = posAtom.x;
		buffer_trafo.getY()[bufIdx] = posAtom.y;
		buffer_trafo.getZ()[bufIdx] = posAtom.z;

	}
}



template<typename REAL>
__global__ void d_rotateForces(
		d_Protein<REAL> const&  protein,
		DOF_6D_Modes<REAL>* dofs,
		unsigned const numDofs
)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned const numAtoms = protein.numAtoms;

	if (idx < protein.numAtoms * numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / protein.numAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;

		Vec3<REAL> ForceAtom(xForce[atomIdx], yForce[atomIdx], zForce[atomIdx]);
		const RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);

		ForceAtom=rotMat*ForceAtom;

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
		DOF_T<REAL>* dofs,
		unsigned const numDOFs,
		Buffer<REAL> buffer_defo,
		Buffer<REAL> buffer_trafo
		)
{
	cudaVerifyKernel((
		d_DOFPos<<<gridSize, blockSize, 0, stream>>> (
		blockSize,
		gridSize,
		stream,
		protein,
		dofs,
		numDOFs,
		buffer_defo,
		buffer_trafo
		))
		);
}



template
 void d_DOFPos<float, DOF_6D_Modes<float>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D_Modes<float>* dofs,
		unsigned const numDOFs,
		Buffer<float> buffer_defo,
		Buffer<float> buffer_trafo
		);

template
 void d_DOFPos<double, DOF_6D_Modes<double>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D_Modes<double>* dofs,
		unsigned const numDOFs,
		Buffer<double> buffer_defo,
		Buffer<double> buffer_trafo
		);

template
 void d_DOFPos<float, DOF_6D_Modes<float>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D_Modes<float>* dofs,
		unsigned const numDOFs,
		Buffer<float> buffer_defo,
		Buffer<float> buffer_trafo
		);

template
 void d_DOFPos<double, DOF_6D_Modes<double>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D_Modes<double>* dofs,
		unsigned const numDOFs,
		Buffer<double> buffer_defo,
		Buffer<double> buffer_trafo
		);

template
 void d_DOFPos<float, DOF_6D<float>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D<float>* dofs,
		unsigned const numDOFs,
		Buffer<float> buffer_defo,
		Buffer<float> buffer_trafo
		);

template
 void d_DOFPos<double, DOF_6D<double>, 0>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D<double>* dofs,
		unsigned const numDOFs,
		Buffer<double> buffer_defo,
		Buffer<double> buffer_trafo
		);

template
 void d_DOFPos<float, DOF_6D<float>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const*  protein,
		DOF_6D<float>* dofs,
		unsigned const numDOFs,
		Buffer<float> buffer_defo,
		Buffer<float> buffer_trafo
		);

template
 void d_DOFPos<double, DOF_6D<double>, 1>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const*  protein,
		DOF_6D<double>* dofs,
		unsigned const numDOFs,
		Buffer<double> buffer_defo,
		Buffer<double> buffer_trafo
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

}  // namespace as
