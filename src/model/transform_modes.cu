#include "transform_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"

namespace as {

/**
 *
 * wouldnt it be better to load consecutive arrays of position into each thread? positions are not consequtive ->guess it takes really long to load
 * or does it average out across each warp?
 */

///



template<typename REAL>
inline __device__ Vec3<REAL> invertDOF	(Vec3<REAL> posAtom,const RotMat<REAL> rotMat){
	const RotMat<REAL> rotMatInv=rotMat.getInv();
	Vec3<REAL> posInv=rotMatInv*dof._6D.pos.inv();
return posInv;
}


template<typename REAL>
void d_deform	(Vec3<REAL>& posAtom,int atomIdx, REAL* dMode, REAL* xModes,REAL* yModes,REAL* zModes, int numModes){
	for(int mode=0;mode<numModes;mode++){
		posAtom.x+=dMode[mode]*xModes[idx*numModes+mode];
		posAtom.y+=dMode[mode]*yModes[idx*numModes+mode];
		posAtom.z+=dMode[mode]*zModes[idx*numModes+mode];
	}
}


template<typename REAL>
void d_translate_rotate	(Vec3<REAL>& posAtom, int atomIdx, Vec3<REAL> const& pos,const RotMat<REAL> rotMat){
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
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFsLig,
		REAL* xRecDefo,
		REAL* yRecDefo,
		REAL* zRecDefo,
		REAL* xRecTrafo,
		REAL* yRecTrafo,
		REAL* yRecTrafo,
		REAL* xLigTrafo,
		REAL* yLigTrafo,
		REAL* zLigTrafo
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

	if (idx < max(numAtomsLig,numAtomsRec)*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;

		const RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);

		if (atomIdx < numDOFsRec){
			Vec3<REAL> posAtomRec(xRec[atomIdx], yRec[atomIdx], zRec[atomIdx]);
			d_deform( posAtomRec,atomIdx, dof.modesRec,  xModesRec, yModesRec, zModesRec, numModesRec);

			xRecDefo[idx] = posAtomRec.x;
			yRecDefo[idx] = posAtomRec.y;
			zRecDefo[idx] = posAtomRec.z;

			Vec3<REAL> posInv=invertDOF(posAtomRec, rotMat);
			d_translate_rotate(posAtomRec,atomIdx, posInv, rotMat.getInv());

			xRecTrafo[idx] = posAtomRec.x;
			xRecTrafo[idx] = posAtomRec.y;
			xRecTrafo[idx] = posAtomRec.z;
		}

		if (atomIdx < numDOFsLig){
			Vec3<REAL> posAtomLig(xLig[atomIdx], yLig[atomIdx], zLig[atomIdx]);
			d_deform( posAtomLig, atomIdx, dof.modesLig,  xModesLig, yModesLig, zModesLig, numModesLig);
			d_translate_rotate(posAtomLig, atomIdx, posAtomLig, rotMat);

			xLigTrafo[idx] = posAtomLig.x;
			xLigTrafo[idx] = posAtomLig.y;
			xLigTrafo[idx] = posAtomLig.z;
		}

	}
}


template<typename REAL>
__global__ void d_rotateForces(
		REAL const* xForce,
		REAL const* yForce,
		REAL const* zForce,
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
		REAL const* xForce,
		REAL const* yForce,
		REAL const* zForce,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numDOFs
)
{
	d_DOFPos<<<gridSize, blockSize, 0, stream>>> (
			xForce,
			yForce,
			zForce,
			dofs,
			numAtoms,
			numDOFs
			);
}


template
void d_rotateForces(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		float const* xForce,
		float const* yForce,
		float const* zForce,
		DOF_6D_Modes<float>* dofs,
		unsigned numAtoms,
		unsigned numDOFs
		);

template
void d_rotateForces(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		double const* xForce,
		double const* yForce,
		double const* zForce,
		DOF_6D_Modes<double>* dofs,
		unsigned numAtoms,
		unsigned numDOFs
		);

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
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		REAL* Trafo,
		REAL* yRecDefo,
		REAL* zRecDefo,
		REAL* xRecTrafo,
		REAL* yRecTrafo,
		REAL* zRecTrafo,
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
				numModesRec,
				numModesLig,
				numDOFs,
				xRecDefo,
				yRecDefo,
				zRecDefo,
				xRecDefo,
				yRecTrafo,
				zRecTrafo,
				xLigTrafo,
				yLigTrafo,
				zLigTrafo
			)
		));
}
template
void d_DOFPos<float>(
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
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		float* xRecDefo,
		float* yRecDefo,
		float* zRecDefo,
		float* xRecTrafo,
		float* yRecTrafo,
		float* zRecTrafo,
		float* xLigTrafo,
		float* yLigTrafo,
		float* zLigTrafo
		);

template
void d_DOFPos<double>(
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
		unsigned numModesRec,
		unsigned numModesLig,
		unsigned numDOFs,
		double* xRecDefo,
		double* yRecDefo,
		double* zRecDefo,
		double* xRecTrafo,
		double* yRecTrafo,
		double* zRecTrafo,
		double* xLigTrafo,
		double* yLigTrafo,
		double* zLigTraf
		);


}  // namespace as
