/*
 * reduction_MB_modes.h
 *
 *  Created on: Nov 22, 2017
 *      Author: glenn
 */

#ifndef REDUCTION_MB_MODES_H_
#define REDUCTION_MB_MODES_H_

#include "Types_6D_MB_Modes.h"
#include "reduction_modes.h"

namespace as {









#ifdef CUDA

template <class T>
void d_reduceMode(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtomsRec,
		const unsigned& numAtomsLig,
		const unsigned& numModesRec,
		const unsigned& numModesLig,
		DOF_6D_MB_Modes<T>* dofs,
		unsigned ligIdx,
		T* xPos, T* yPos, T* zPos,
		T *xModesLig,T *yModesLig,T *zModesLig,
		T *xModesRec,T *yModesRec,T *zModesRec,
		T *d_fxRec, T *d_fyRec, T *d_fzRec,
		T *d_fxLig, T *d_fyLig, T *d_fzLig,
		T *d_E,
		T *g_odata,
		const cudaStream_t& stream);

/* wrapper function to call the device kernel for the reduction */
template<typename REAL>
void deviceReduce(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned& receptorSize,
		const unsigned& ligandSize,
		const unsigned& numModesRec,
		const unsigned& numModesLig,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned ligIdx,
		REAL *xModesRec,REAL *yModesRec,REAL *zModesRec,
		REAL *xModesLig,REAL *yModesLig,REAL *zModesLig,
		REAL* xPos, REAL* yPos, REAL* zPos,
		REAL *d_fxRec, REAL *d_fyRec, REAL *d_fzRec,
		REAL *d_fxLig, REAL *d_fyLig, REAL *d_fzLig,
		REAL *d_E,
		REAL *d_out,
		const cudaStream_t& stream)
{
	/* we need at least twice the block size number of threads */
	//const unsigned threads = (max(ligandSize, receptorSize) < blockSize*2) ? nextPow2((max(ligandSize, receptorSize) + 1)/ 2) : blockSize;
	unsigned size = max(ligandSize, receptorSize);
	const unsigned threads = (size < blockSize*2) ? nextPow2((size + 1)/ 2) : blockSize;
	/* each structure is reduced by one thread block */
	const unsigned blocks = numDOFs;
	d_reduceMode(threads, blocks, receptorSize, ligandSize, numModesRec, numModesLig,
			dofs, ligIdx,
			xPos, yPos, zPos,
			xModesRec, yModesRec, zModesRec,
			xModesLig, yModesLig, zModesLig,
			d_fxRec, d_fyRec, d_fzRec,
			d_fxLig, d_fyLig, d_fzLig,
			d_E,
			d_out,
			stream);
}

/* remaining reduce part that is performed on the host */
template<typename REAL>
void h_finalReduce(
			const unsigned& numDOFs,
			DOF_6D_MB_Modes<REAL>* dofs,
			unsigned int ligIdx,
			REAL const* modeForceConstantRec,
			REAL const* modeForceConstantLig,
			const unsigned int& numModesRec,
			const unsigned int& numModesLig,
			const REAL* deviceOut,
			Result_6D_MB_Modes<REAL>* enGrads)
{
	unsigned const dofSize = 13 + numModesRec + numModesLig;
	for (unsigned i = 0; i < numDOFs; ++i)
	{
		auto &enGrad = enGrads[i];
		enGrad._6D[ligIdx].pos.x = deviceOut[i*dofSize + 0];
		enGrad._6D[ligIdx].pos.y = deviceOut[i*dofSize + 1];
		enGrad._6D[ligIdx].pos.z = deviceOut[i*dofSize + 2];

		for(unsigned j = 0; j < 3; ++j) {
			REAL magn2 = enGrad._6D[ligIdx].pos.x*enGrad._6D[ligIdx].pos.x
					+ enGrad._6D[ligIdx].pos.y*enGrad._6D[ligIdx].pos.y
					+ enGrad._6D[ligIdx].pos.z*enGrad._6D[ligIdx].pos.z;

			if(magn2 > static_cast<REAL>(ForceLim)) {
				enGrad._6D[ligIdx].pos.x *= 0.01;
				enGrad._6D[ligIdx].pos.y *= 0.01;
				enGrad._6D[ligIdx].pos.z *= 0.01;
			}
		}

		enGrad.E = deviceOut[i*dofSize + 3];

		Torque<REAL> torque;
		torque.mat[0][0] = deviceOut[i*dofSize + 4 ];
		torque.mat[0][1] = deviceOut[i*dofSize + 5 ];
		torque.mat[0][2] = deviceOut[i*dofSize + 6 ];
		torque.mat[1][0] = deviceOut[i*dofSize + 7 ];
		torque.mat[1][1] = deviceOut[i*dofSize + 8 ];
		torque.mat[1][2] = deviceOut[i*dofSize + 9 ];
		torque.mat[2][0] = deviceOut[i*dofSize + 10];
		torque.mat[2][1] = deviceOut[i*dofSize + 11];
		torque.mat[2][2] = deviceOut[i*dofSize + 12];

		for(int mode=0; mode < numModesLig; mode++){
			enGrad.modesLig[ligIdx][mode]=deviceOut[i*dofSize + 13 +mode];
		}

		for(int mode=0; mode < numModesRec; mode++){
			enGrad.modesRec[mode]=deviceOut[i*dofSize + 13 + numModesLig + mode];
		}

		correctModeForce(
			modeForceConstantRec,
			numModesRec,
			enGrad.modesRec
			);

		correctModeForce(
			modeForceConstantLig,
			numModesLig,
			enGrad.modesLig[ligIdx]
			);


		const auto &dof = dofs[i];
		const TorqueMat<REAL> torqueMat = euler2torquemat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);
		Vec3<REAL> result = torqueMat.rotateReduce(torque);

		enGrad._6D[ligIdx].ang.x = result.x;
		enGrad._6D[ligIdx].ang.y = result.y;
		enGrad._6D[ligIdx].ang.z = result.z;
	}
}
#endif // CUDA





}//end namespace as
#endif /* REDUCTION_MODES_H_ */
