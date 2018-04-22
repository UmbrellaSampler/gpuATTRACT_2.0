#include "transform_modes.h"

#include "Vec3.h"
#include "RotMat.h"
#include "matrixFunctions.h"
#include "macros.h"
namespace as {




template<typename REAL, bool PROTEINTYPE, bool MODES>
__global__ void d_DOFPos(
		REAL const* xPos,
		REAL const* yPos,
		REAL const* zPos,
		REAL const* xModes,
		REAL const* yModes,
		REAL const* zModes,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numModes,
		unsigned numDOFs,
		REAL* xDefo,
		REAL* yDefo,
		REAL* zDefo,
		REAL* xTrafo,
		REAL* yTrafo,
		REAL* zTrafo
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < numAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;

		Vec3<REAL> 		posAtom(xPos[atomIdx], yPos[atomIdx], zPos[atomIdx]);
		RotMat<REAL> 	rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
		Vec3<REAL> 		translation(dof._6D.pos.x, dof._6D.pos.y, dof._6D.pos.z);
		if ( idx % numAtoms == 0){

			//printf("%f %f %f %f  \n modes r %f %f %f %f %f  %d\n modes l %f %f %f %f %f %d\n",xModes[atomIdx*numModes], dof._6D.ang.y, dof._6D.ang.z,xPos[atomIdx], dof.modesRec[0],dof.modesRec[1],dof.modesRec[2],dof.modesRec[3],dof.modesRec[4], DOFidx, dof.modesLig[0],dof.modesLig[1],dof.modesLig[2],dof.modesLig[3],dof.modesLig[4], DOFidx);
		}
		if(MODES == true){
			if(PROTEINTYPE == false){
				rotMat = rotMat.getInv();
				translation = rotMat * translation.inv();
				for(int mode=0; mode < 1; mode++){
					posAtom.x += dof.modesRec[mode] * xModes[atomIdx*numModes+mode];
					posAtom.y += dof.modesRec[mode] * yModes[atomIdx*numModes+mode];
					posAtom.z += dof.modesRec[mode] * zModes[atomIdx*numModes+mode];
//					posAtom.x += dof.modesRec[0];
//					posAtom.y += dof.modesRec[0];
//					posAtom.z += dof.modesRec[0];
					//printf("modes 2 rec %f %d %d %d\n",dof.modesRec[0],PROTEINTYPE,DOFidx,atomIdx);
//					posAtom.x += dof.modesRec[mode] ;
//					posAtom.y += dof.modesRec[mode] ;
//					posAtom.z += dof.modesRec[mode] ;
				}
			}
			else {
				for(int mode=0; mode < numModes; mode++){
//					posAtom.x += dof.modesLig[mode] * xModes[atomIdx*numModes+mode];
//					posAtom.y += dof.modesLig[mode] * yModes[atomIdx*numModes+mode];
//					posAtom.z += dof.modesLig[mode] * zModes[atomIdx*numModes+mode];

//					posAtom.x += dof.modesLig[mode] * xModes[atomIdx*numModes+mode];
//					posAtom.y += dof.modesLig[mode] * yModes[atomIdx*numModes+mode];
//					posAtom.z += dof.modesLig[mode] * zModes[atomIdx*numModes+mode];
				}
			}

			xDefo[idx] = posAtom.x;
			yDefo[idx] = posAtom.y;
			zDefo[idx] = posAtom.z;
		}
		posAtom = rotMat * posAtom;
		posAtom += translation;

		xTrafo[idx] = posAtom.x;
		yTrafo[idx] = posAtom.y;
		zTrafo[idx] = posAtom.z;
	}
}



template<typename REAL, bool MODES>
__global__ void d_DOFPos_rec(
		REAL const* xPos,
		REAL const* yPos,
		REAL const* zPos,
		REAL const* xModes,
		REAL const* yModes,
		REAL const* zModes,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numModes,
		unsigned numDOFs,
		REAL* xDefo,
		REAL* yDefo,
		REAL* zDefo,
		REAL* xTrafo,
		REAL* yTrafo,
		REAL* zTrafo
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < numAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;

		Vec3<REAL> 		posAtom(xPos[atomIdx], yPos[atomIdx], zPos[atomIdx]);
		RotMat<REAL> 	rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
		Vec3<REAL> 		translation(dof._6D.pos.x, dof._6D.pos.y, dof._6D.pos.z);



			rotMat = rotMat.getInv();
			translation = rotMat * translation.inv();
			for(int mode=0; mode < numModes; mode++){
				posAtom.x += dof.modesRec[mode] * xModes[atomIdx*numModes+mode];
				posAtom.y += dof.modesRec[mode] * yModes[atomIdx*numModes+mode];
				posAtom.z += dof.modesRec[mode] * zModes[atomIdx*numModes+mode];
			}


		xDefo[idx] = posAtom.x;
		yDefo[idx] = posAtom.y;
		zDefo[idx] = posAtom.z;

		posAtom = rotMat * posAtom;
		posAtom += translation;

		xTrafo[idx] = posAtom.x;
		yTrafo[idx] = posAtom.y;
		zTrafo[idx] = posAtom.z;
	}
}

template<typename REAL, bool MODES>
__global__ void d_DOFPos_lig(
		REAL const* xPos,
		REAL const* yPos,
		REAL const* zPos,
		REAL const* xModes,
		REAL const* yModes,
		REAL const* zModes,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numAtoms,
		unsigned numModes,
		unsigned numDOFs,
		REAL* xDefo,
		REAL* yDefo,
		REAL* zDefo,
		REAL* xTrafo,
		REAL* yTrafo,
		REAL* zTrafo
		)
{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < numAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;

		Vec3<REAL> 		posAtom(xPos[atomIdx], yPos[atomIdx], zPos[atomIdx]);
		RotMat<REAL> 	rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
		Vec3<REAL> 		translation(dof._6D.pos.x, dof._6D.pos.y, dof._6D.pos.z);


		for(int mode=0; mode < numModes; mode++){
			posAtom.x += dof.modesLig[mode] * xModes[atomIdx*numModes+mode];
			posAtom.y += dof.modesLig[mode] * yModes[atomIdx*numModes+mode];
			posAtom.z += dof.modesLig[mode] * zModes[atomIdx*numModes+mode];
		}
		xDefo[idx] = posAtom.x;
		yDefo[idx] = posAtom.y;
		zDefo[idx] = posAtom.z;

		posAtom = rotMat * posAtom;
		posAtom += translation;

		xTrafo[idx] = posAtom.x;
		yTrafo[idx] = posAtom.y;
		zTrafo[idx] = posAtom.z;
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

		Vec3<REAL> ForceAtom(xForce[idx], yForce[idx], zForce[idx]);
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




template<typename REAL, bool PROTEINTYPE, bool MODES>
void d_DOFPos(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL>* protein,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numDOFs,
		REAL* xDefo,
		REAL* yDefo,
		REAL* zDefo,
		REAL* xTrafo,
		REAL* yTrafo,
		REAL* zTrafo
		){
	cudaVerifyKernel((
			d_DOFPos<REAL,PROTEINTYPE, MODES><<<gridSize, blockSize, 0, stream>>> (
				protein->xPos,
				protein->yPos,
				protein->zPos,
				protein->xModes,
				protein->yModes,
				protein->zModes,
				dofs,
				protein->numAtoms,
				protein->numModes,
				numDOFs,
				xDefo,
				yDefo,
				zDefo,
				xTrafo,
				yTrafo,
				zTrafo
			))
		);
}

template
void d_DOFPos<float, false, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos<double, false, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template
void d_DOFPos<float, true, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos<double, true, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template
void d_DOFPos<float, false, true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos<double, false, true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template
void d_DOFPos<float, true, true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos<double, true, true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
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

///////////////////////////////////////////////

template
void d_DOFPos_lig<float, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos_lig<double,  false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template
void d_DOFPos_lig<float, true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos_lig<double, true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template
void d_DOFPos_rec<float,  true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos_rec<double,  true >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template
void d_DOFPos_rec<float, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		unsigned numDOFs,
		float* xDefo,
		float* yDefo,
		float* zDefo,
		float* xTrafo,
		float* yTrafo,
		float* zTrafo
		);

template
void d_DOFPos_rec<double, false >(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		unsigned numDOFs,
		double* xDefo,
		double* yDefo,
		double* zDefo,
		double* xTrafo,
		double* yTrafo,
		double* zTrafo
		);

template<typename REAL,  bool MODES>
void d_DOFPos_rec(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL>* protein,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numDOFs,
		REAL* xDefo,
		REAL* yDefo,
		REAL* zDefo,
		REAL* xTrafo,
		REAL* yTrafo,
		REAL* zTrafo
		){
	cudaVerifyKernel((
			d_DOFPos_rec<REAL, MODES><<<gridSize, blockSize, 0, stream>>> (
				protein->xPos,
				protein->yPos,
				protein->zPos,
				protein->xModes,
				protein->yModes,
				protein->zModes,
				dofs,
				protein->numAtoms,
				protein->numModes,
				numDOFs,
				xDefo,
				yDefo,
				zDefo,
				xTrafo,
				yTrafo,
				zTrafo
			))
		);
}
template<typename REAL, bool MODES>
void d_DOFPos_lig(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL>* protein,
		DOF_6D_Modes<REAL>* dofs,
		unsigned numDOFs,
		REAL* xDefo,
		REAL* yDefo,
		REAL* zDefo,
		REAL* xTrafo,
		REAL* yTrafo,
		REAL* zTrafo
		){
	cudaVerifyKernel((
			d_DOFPos_lig<REAL, MODES><<<gridSize, blockSize, 0, stream>>> (
				protein->xPos,
				protein->yPos,
				protein->zPos,
				protein->xModes,
				protein->yModes,
				protein->zModes,
				dofs,
				protein->numAtoms,
				protein->numModes,
				numDOFs,
				xDefo,
				yDefo,
				zDefo,
				xTrafo,
				yTrafo,
				zTrafo
			))
		);
}

}  // namespace as
