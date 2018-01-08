/*
 * transform_MB_modes.h
 *
 *  Created on: Dec 5, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TRANSFORM_MB_MODES_H_
#define SRC_MODEL_TRANSFORM_MB_MODES_H_


#include "transform_modes.h"
#include "Types_6D_MB_Modes.h"


 namespace as{


 template<typename REAL>
 const DOF_6D_MB_Modes<REAL> invertDOF( DOF_6D_MB_Modes<REAL> ligandDOF, const unsigned int numLigands)
 {
 	DOF_6D_MB_Modes<REAL> invertedDOF;
 	for (unsigned int lig = 0; lig < numLigands; lig++){
		Vec3<REAL> ang(0.0);
		invertedDOF._6D[lig].ang=ligandDOF._6D[lig].ang.inv();
		invertedDOF._6D[lig].pos=ligandDOF._6D[lig].pos.inv();
		const RotMat<REAL> rotMat=euler2rotmat(ligandDOF._6D[lig].ang.x, ligandDOF._6D[lig].ang.y, ligandDOF._6D[lig].ang.z).getInv();
		invertedDOF._6D[lig].pos=rotMat*invertedDOF._6D[lig].pos;

		std::copy(ligandDOF.modesLig[lig], ligandDOF.modesLig[lig]+MODES_MAX_NUMBER,
			invertedDOF.modesLig[lig]);
 	}
 	std::copy( ligandDOF.modesRec, ligandDOF.modesRec+MODES_MAX_NUMBER,
			invertedDOF.modesRec);
 return invertedDOF;
 }




#ifdef CUDA

 template<typename REAL>
 void d_DOFPos(
 		unsigned blockSize,
 		unsigned gridSize,
 		const cudaStream_t &stream,
 		unsigned numLigands,
 		unsigned ligIdx,
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
 		DOF_6D_MB_Modes<REAL>* dofs,
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
 		REAL** xLigTrafo,
 		REAL** yLigTrafo,
 		REAL** zLigTrafo
 		);

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
);

#endif


}

#endif /* TRANSFORM_MODES_H_ */
