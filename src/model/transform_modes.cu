#include "transform_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"
namespace as {




template<typename REAL>
inline __device__ Vec3<REAL> invertDOF	(Vec3<REAL> posAtom,const RotMat<REAL> rotMat, Vec3<REAL> pos){
	const RotMat<REAL> rotMatInv = rotMat.getInv();
	Vec3<REAL> posInv = rotMatInv * pos.inv();
return posInv;
}


template<typename REAL>
inline __device__ void d_deform(Vec3<REAL>& posAtom,int atomIdx, REAL* dMode, REAL* xModes,REAL* yModes,REAL* zModes, int numModes){
	for(int mode=0; mode < numModes; mode++){
		posAtom.x += dMode[mode] * xModes[atomIdx*numModes+mode];
		posAtom.y += dMode[mode] * yModes[atomIdx*numModes+mode];
		posAtom.z += dMode[mode] * zModes[atomIdx*numModes+mode];
	}
}

template<typename REAL>
inline __device__ void d_translate_rotate	(Vec3<REAL>& posAtom, int atomIdx, Vec3<REAL> const& pos,const RotMat<REAL> rotMat){
	posAtom = rotMat*posAtom;
	posAtom += pos;
}


template<typename REAL>
__global__ void d_DOFPos(
		REAL const* xRec,
		REAL const* yRec,
		REAL const* zRec,
		REAL const* xLig,
		REAL const* yLig,
		REAL const* zLig,
		REAL const* xModesRec,
		REAL const* yModesRec,
		REAL const* zModesRec,
		REAL const* xModesLig,
		REAL const* yModesLig,
		REAL const* zModesLig,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtomsRec,
		unsigned numAtomsLig,
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		REAL* xRecDefo,
		REAL* yRecDefo,
		REAL* zRecDefo,
		REAL* xRecTrafo,
		REAL* yRecTrafo,
		REAL* zRecTrafo,
		REAL* xLigDefo,
		REAL* yLigDefo,
		REAL* zLigDefo,
		REAL* xLigTrafo,
		REAL* yLigTrafo,
		REAL* zLigTrafo
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned int maxNumAtoms = max(numAtomsRec, numAtomsLig);


	if (idx < maxNumAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / maxNumAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % maxNumAtoms;

		const RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);

		if (atomIdx < numAtomsRec ){
			int bufIdx = numAtomsRec * DOFidx + atomIdx;
			Vec3<REAL> posAtomRec(xRec[atomIdx], yRec[atomIdx], zRec[atomIdx]);

			for(int mode=0; mode < numModesRec; mode++){
				posAtomRec.x += dof.modesRec[mode] * xModesRec[atomIdx*numModesRec+mode];
				posAtomRec.y += dof.modesRec[mode] * yModesRec[atomIdx*numModesRec+mode];
				posAtomRec.z += dof.modesRec[mode] * zModesRec[atomIdx*numModesRec+mode];
			}

			xRecDefo[bufIdx] = posAtomRec.x;
			yRecDefo[bufIdx] = posAtomRec.y;
			zRecDefo[bufIdx] = posAtomRec.z;

			const RotMat<REAL> rotMatInv = rotMat.getInv();
			Vec3<REAL> posInv = rotMatInv * dof._6D.pos.inv();
			posAtomRec = rotMatInv*posAtomRec;
			posAtomRec += posInv;


			xRecTrafo[bufIdx] = posAtomRec.x;
			yRecTrafo[bufIdx] = posAtomRec.y;
			zRecTrafo[bufIdx] = posAtomRec.z;
		}

		if (atomIdx < numAtomsLig && idx < numAtomsRec*numDOFs){
			int bufIdx = numAtomsLig * DOFidx + atomIdx;

			Vec3<REAL> posAtomLig(xLig[atomIdx], yLig[atomIdx], zLig[atomIdx]);


			for(int mode=0; mode < numModesLig; mode++){
				posAtomLig.x += dof.modesLig[mode] * xModesLig[atomIdx*numModesLig+mode];
				posAtomLig.y += dof.modesLig[mode] * yModesLig[atomIdx*numModesLig+mode];
				posAtomLig.z += dof.modesLig[mode] * zModesLig[atomIdx*numModesLig+mode];
			}

			xLigDefo[bufIdx] = posAtomLig.x;
			yLigDefo[bufIdx] = posAtomLig.y;
			zLigDefo[bufIdx] = posAtomLig.z;

			posAtomLig = rotMat*posAtomLig;
			posAtomLig += dof._6D.pos;

			xLigTrafo[bufIdx] = posAtomLig.x;
			yLigTrafo[bufIdx] = posAtomLig.y;
			zLigTrafo[bufIdx] = posAtomLig.z;
		}
	}
}


template<typename REAL>
__global__ void d_rotateForces(
		REAL* xForce,
		REAL* yForce,
		REAL* zForce,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numDOFs
)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < numAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
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




template<typename REAL>
void d_DOFPos(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		REAL const* xRec,
		REAL const* yRec,
		REAL const* zRec,
		REAL const* xLig,
		REAL const* yLig,
		REAL const* zLig,
		REAL const* xModesRec,
		REAL const* yModesRec,
		REAL const* zModesRec,
		REAL const* xModesLig,
		REAL const* yModesLig,
		REAL const* zModesLig,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtomsRec,
		unsigned numAtomsLig,
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		REAL* xRecDefo,
		REAL* yRecDefo,
		REAL* zRecDefo,
		REAL* xRecTrafo,
		REAL* yRecTrafo,
		REAL* zRecTrafo,
		REAL* xLigDefo,
		REAL* yLigDefo,
		REAL* zLigDefo,
		REAL* xLigTrafo,
		REAL* yLigTrafo,
		REAL* zLigTrafo)
{
	cudaVerifyKernel((
			d_DOFPos<<<gridSize, blockSize, 0, stream>>> (
			xRec,
			yRec,
			zRec,
			xLig,
			yLig,
			zLig,
			xModesRec,
			yModesRec,
			zModesRec,
			xModesLig,
			yModesLig,
			zModesLig,
			dofs,
			numAtomsRec,
			numAtomsLig,
			numModesRec,
			numModesLig,
			numDOFs,
			xRecDefo,
			yRecDefo,
			zRecDefo,
			xRecTrafo,
			yRecTrafo,
			zRecTrafo,
			xLigDefo,
			yLigDefo,
			zLigDefo,
			xLigTrafo,
			yLigTrafo,
			zLigTrafo
			))
		);
}

template
void d_DOFPos<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		float const* xRec,
		float const* yRec,
		float const* zRec,
		float const* xLig,
		float const* yLig,
		float const* zLig,
		float const* xModesRec,
		float const* yModesRec,
		float const* zModesRec,
		float const* xModesLig,
		float const* yModesLig,
		float const* zModesLig,
		DOF_6D_Modes<float>* dofs,
		unsigned numAtomsRec,
		unsigned numAtomsLig,
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		float* xRecDefo,
		float* yRecDefo,
		float* zRecDefo,
		float* xRecTrafo,
		float* yRecTrafo,
		float* zRecTrafo,
		float* xLigDefo,
		float* yLigDefo,
		float* zLigDefo,
		float* xLigTrafo,
		float* yLigTrafo,
		float* zLigTrafo);

template
void d_DOFPos<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		double const* xRec,
		double const* yRec,
		double const* zRec,
		double const* xLig,
		double const* yLig,
		double const* zLig,
		double const* xModesRec,
		double const* yModesRec,
		double const* zModesRec,
		double const* xModesLig,
		double const* yModesLig,
		double const* zModesLig,
		DOF_6D_Modes<double>* dofs,
		unsigned numAtomsRec,
		unsigned numAtomsLig,
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		double* xRecDefo,
		double* yRecDefo,
		double* zRecDefo,
		double* xRecTrafo,
		double* yRecTrafo,
		double* zRecTrafo,
		double* xLigDefo,
		double* yLigDefo,
		double* zLigDefo,
		double* xLigTrafo,
		double* yLigTrafo,
		double* zLigTrafo);



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
