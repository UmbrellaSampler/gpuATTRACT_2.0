#include "transform_MB_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"
namespace as {


/**
 * This function is ment for multi body transformation and deformation. It calculates all necessary positions for one ligand in relation to all other ligands in the system
 * and the receptor. Assuming that the partner of this ligand is located in the origin.
 *
 */
template<typename REAL>
__global__ void d_DOFPos(
		unsigned numLigands,
		unsigned ligIdx,
		REAL const* xLig,
		REAL const* yLig,
		REAL const* zLig,
		REAL const* xModesLig,
		REAL const* yModesLig,
		REAL const* zModesLig,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned numAtomsLig,
		unsigned numModesLig,
		unsigned numDOFs,
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

	/* load DOF from global memory */
	unsigned DOFidx = idx / numAtomsLig;
	auto dof = dofs[DOFidx];
	unsigned atomIdx = idx % numAtomsLig;

		if (idx < numAtomsLig * numDOFs) {
			const RotMat<REAL> rotMat = euler2rotmat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);
			unsigned bufSize = numDOFs * numAtomsLig *numLigands;
			int bufIdx = numAtomsLig * DOFidx + atomIdx;

			Vec3<REAL> posAtomLig(xLig[atomIdx], yLig[atomIdx], zLig[atomIdx]);


			for(int mode=0; mode < numModesLig; mode++){
				posAtomLig.x += dof.modesLig[ligIdx][mode] * xModesLig[atomIdx*numModesLig+mode];
				posAtomLig.y += dof.modesLig[ligIdx][mode] * yModesLig[atomIdx*numModesLig+mode];
				posAtomLig.z += dof.modesLig[ligIdx][mode] * zModesLig[atomIdx*numModesLig+mode];
			}

			xLigDefo[bufIdx] = posAtomLig.x;
			yLigDefo[bufIdx] = posAtomLig.y;
			zLigDefo[bufIdx] = posAtomLig.z;


			for( int lig = 0; lig < numLigands; lig++){
				unsigned bufIdxMB = lig * bufSize + bufIdx;
				if( lig != ligIdx){
					//get the invers of the rotationmatrix of each ligand
					const RotMat<REAL> rotMatInv = euler2rotmat(dof._6D[lig].ang.x, dof._6D[lig].ang.y, dof._6D[lig].ang.z).getInv();
					//get the relative positon of ligand[ligIdx] to ligang[lig]
					Vec3<REAL> tRel =  dof._6D[ligIdx].pos - dof._6D[lig].pos;
					//rotate tRel into the system of ligand[lig]
					tRel = rotMatInv * tRel;
					//rotate each position of ligand[ligIdx] into the system of ligand[lig] and than rotate by the angle of ligand[ligIdx]
					posAtomLig = rotMatInv * posAtomLig;
					posAtomLig = rotMat *  posAtomLig;
					// add the relative translation with is now in the coordinate system of ligand[lig]
					posAtomLig += tRel;

					xLigTrafo[bufIdxMB] = posAtomLig.x;
					yLigTrafo[bufIdxMB] = posAtomLig.y;
					zLigTrafo[bufIdxMB] = posAtomLig.z;
				}
				//if lig == ligIdx calculate the orientation of ligand[ligIdx] relative to the receptor
				else if( lig == ligIdx){
					posAtomLig = rotMat*posAtomLig;
					posAtomLig += dof._6D[lig].pos;

					xLigTrafo[bufIdxMB] = posAtomLig.x;
					yLigTrafo[bufIdxMB] = posAtomLig.y;
					zLigTrafo[bufIdxMB] = posAtomLig.z;
				}
			}
		}

}


template<typename REAL>
__global__ void d_rotateForces(
		REAL* xForce,
		REAL* yForce,
		REAL* zForce,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		unsigned ligIdx
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
		const RotMat<REAL> rotMat = euler2rotmat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);

		ForceAtom=rotMat*ForceAtom;

		xForce[idx] = ForceAtom.x;
		yForce[idx] = ForceAtom.y;
		zForce[idx] = ForceAtom.z;
	}
}

template<typename REAL>
__global__ void d_DOFPos_receptor(
		unsigned numLigands,
		REAL const* xRec,
		REAL const* yRec,
		REAL const* zRec,
		REAL const* xModesRec,
		REAL const* yModesRec,
		REAL const* zModesRec,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned numAtomsRec,
		unsigned numModesRec,
		unsigned numDOFs,
		REAL* xRecDefo,
		REAL* yRecDefo,
		REAL* zRecDefo,
		REAL* xRecTrafo,
		REAL* yRecTrafo,
		REAL* zRecTrafo
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;



	if (idx < numAtomsRec * numDOFs) {
		unsigned DOFidx = idx / numAtomsRec;
		unsigned atomIdx = idx % numAtomsRec;
		unsigned bufSize = numAtomsRec * numLigands * numDOFs;
		int bufIdx = numAtomsRec * DOFidx + atomIdx;
		auto dof = dofs[DOFidx];


		Vec3<REAL> posAtomRec(xRec[bufIdx], yRec[bufIdx], zRec[bufIdx]);

		for(int mode=0; mode < numModesRec; mode++){
			posAtomRec.x += dof.modesRec[mode] * xModesRec[atomIdx*numModesRec+mode];
			posAtomRec.y += dof.modesRec[mode] * yModesRec[atomIdx*numModesRec+mode];
			posAtomRec.z += dof.modesRec[mode] * zModesRec[atomIdx*numModesRec+mode];
		}

		xRecDefo[bufIdx] = posAtomRec.x;
		yRecDefo[bufIdx] = posAtomRec.y;
		zRecDefo[bufIdx] = posAtomRec.z;
		for(unsigned ligIdx = 0; ligIdx < numLigands; ligIdx++){
			unsigned bufIdxMB = ligIdx * bufSize + bufIdx;
			const RotMat<REAL> rotMat = euler2rotmat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);
			const RotMat<REAL> rotMatInv = rotMat.getInv();
			Vec3<REAL> posInv = rotMatInv * dof._6D[ligIdx].pos.inv();

			posAtomRec = rotMatInv*posAtomRec;
			posAtomRec += posInv;

			xRecTrafo[bufIdxMB] = posAtomRec.x;
			yRecTrafo[bufIdxMB] = posAtomRec.y;
			zRecTrafo[bufIdxMB] = posAtomRec.z;
		}
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
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		unsigned ligIdx
)
{
	d_rotateForces<<<gridSize, blockSize, 0, stream>>> (
			xForce,
			yForce,
			zForce,
			dofs,
			numAtoms,
			numDOFs,
			ligIdx
			);
}




template<typename REAL>
void d_DOFPos(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned numLigands,
		unsigned ligIdx,
		REAL const* xLig,
		REAL const* yLig,
		REAL const* zLig,
		REAL const* xModesLig,
		REAL const* yModesLig,
		REAL const* zModesLig,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned numAtomsLig,
		unsigned numModesLig,
		unsigned numDOFs,
		REAL* xLigDefo,
		REAL* yLigDefo,
		REAL* zLigDefo,
		REAL* xLigTrafo,
		REAL* yLigTrafo,
		REAL* zLigTrafo
		)
{
	cudaVerifyKernel((
			d_DOFPos<<<gridSize, blockSize, 0, stream>>> (
				numLigands,
				ligIdx,
				xLig,
				yLig,
				zLig,
				xModesLig,
				yModesLig,
				zModesLig,
				dofs,
				numAtomsLig,
				numModesLig,
				numDOFs,
				xLigDefo,
				yLigDefo,
				zLigDefo,
				xLigTrafo,
				yLigTrafo,
				zLigTrafo
				))
		);
}

template<typename REAL>
void d_DOFPos_receptor(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned numLigands,
		REAL const* xRec,
		REAL const* yRec,
		REAL const* zRec,
		REAL const* xModesRec,
		REAL const* yModesRec,
		REAL const* zModesRec,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned numAtomsRec,
		unsigned numModesRec,
		unsigned numDOFs,
		REAL* xRecDefo,
		REAL* yRecDefo,
		REAL* zRecDefo,
		REAL* xRecTrafo,
		REAL* yRecTrafo,
		REAL* zRecTrafo
		)
{
	cudaVerifyKernel((
			d_DOFPos_receptor<<<gridSize, blockSize, 0, stream>>> (
				numLigands,
				xRec,
				yRec,
				zRec,
				xModesRec,
				yModesRec,
				zModesRec,
				dofs,
				numAtomsRec,
				numModesRec,
				numDOFs,
				xRecDefo,
				yRecDefo,
				zRecDefo,
				xRecTrafo,
				yRecTrafo,
				zRecTrafo
))
		);
}

template
void d_DOFPos<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned numLigands,
		unsigned ligIdx,
		float const* xLig,
		float const* yLig,
		float const* zLig,
		float const* xModesLig,
		float const* yModesLig,
		float const* zModesLig,
		DOF_6D_MB_Modes<float>* dofs,
		unsigned numAtomsLig,
		unsigned numModesLig,
		unsigned numDOFs,
		float* xLigDefo,
		float* yLigDefo,
		float* zLigDefo,
		float* xLigTrafo,
		float* yLigTrafo,
		float* zLigTrafo
		);

template
void d_DOFPos<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned numLigands,
		unsigned ligIdx,
		double const* xLig,
		double const* yLig,
		double const* zLig,
		double const* xModesLig,
		double const* yModesLig,
		double const* zModesLig,
		DOF_6D_MB_Modes<double>* dofs,
		unsigned numAtomsLig,
		unsigned numModesLig,
		unsigned numDOFs,
		double* xLigDefo,
		double* yLigDefo,
		double* zLigDefo,
		double* xLigTrafo,
		double* yLigTrafo,
		double* zLigTrafo
		);

template
void d_DOFPos_receptor<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned numLigands,
		float const* xRec,
		float const* yRec,
		float const* zRec,
		float const* xModesRec,
		float const* yModesRec,
		float const* zModesRec,
		DOF_6D_MB_Modes<float>* dofs,
		unsigned numAtomsRec,
		unsigned numModesRec,
		unsigned numDOFs,
		float* xRecDefo,
		float* yRecDefo,
		float* zRecDefo,
		float* xRecTrafo,
		float* yRecTrafo,
		float* zRecTrafo
		);

template
void d_DOFPos_receptor<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		unsigned numLigands,
		double const* xRec,
		double const* yRec,
		double const* zRec,
		double const* xModesRec,
		double const* yModesRec,
		double const* zModesRec,
		DOF_6D_MB_Modes<double>* dofs,
		unsigned numAtomsRec,
		unsigned numModesRec,
		unsigned numDOFs,
		double* xRecDefo,
		double* yRecDefo,
		double* zRecDefo,
		double* xRecTrafo,
		double* yRecTrafo,
		double* zRecTrafo
		);

template
void d_rotateForces<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		float* xForce,
		float* yForce,
		float* zForce,
		DOF_6D_MB_Modes<float>* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		unsigned ligIdx);

template
void d_rotateForces<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		double* xForce,
		double* yForce,
		double* zForce,
		DOF_6D_MB_Modes<double>* dofs,
		unsigned numAtoms,
		unsigned numDOFs,
		unsigned ligIdx);

}  // namespace as
