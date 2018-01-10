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
 const DOF_6D_MB_Modes<REAL> invertDOF( DOF_6D_MB_Modes<REAL> ligandDOF, unsigned int numLigands)
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


 template<typename REAL>
  void deform(
  		REAL const* x,
  		REAL const* y,
  		REAL const* z,
  		unsigned const& numAtoms,
  		unsigned const& numModes,
  		REAL const* dlig,
  		REAL const* xModes,
  		REAL const* yModes,
  		REAL const* zModes,
  		REAL* xDeformed,
  		REAL* yDeformed,
  		REAL* zDeformed
  		)
 {
	for (unsigned i = 0; i < numAtoms; ++i) {
		Vec3<REAL> posAtom(x[i], y[i], z[i]);
		for(int mode=0;mode<numModes;mode++){
			posAtom.x+=dlig[mode]*xModes[i*numModes+mode];
			posAtom.y+=dlig[mode]*yModes[i*numModes+mode];
			posAtom.z+=dlig[mode]*zModes[i*numModes+mode];
		}
		xDeformed[i]=posAtom.x;
		yDeformed[i]=posAtom.y;
		zDeformed[i]=posAtom.z;
	}




  }
 template<typename REAL>
 void rotate_translate(
 		REAL const* x,
 		REAL const* y,
 		REAL const* z,
 		Vec3<REAL> const& displacementRec,
		Vec3<REAL> const& angRec,
 		Vec3<REAL> const& displacementLig,
 		Vec3<REAL> const& angLig,
 		unsigned const& numAtoms,
 		REAL* xTr,
 		REAL* yTr,
 		REAL* zTr)
 {
 	for (unsigned i = 0; i < numAtoms; ++i) {
		Vec3<REAL> posAtom(x[i], y[i], z[i]);

		const RotMat<REAL> rotMat = euler2rotmat(angLig.x, angLig.y, angLig.z);
		const RotMat<REAL> rotMatInv = euler2rotmat(angRec.x, angRec.y, angRec.z).getInv();
		//get the relative positon of ligand[ligIdx] to ligang[lig]
		Vec3<REAL> tRel =  displacementLig - displacementRec;
		//rotate tRel into the system of ligand[lig]
		tRel = rotMatInv * tRel;
		//rotate each position of ligand[ligIdx] into the system of ligand[lig] and than rotate by the angle of ligand[ligIdx]
		posAtom = rotMatInv * posAtom;
		posAtom = rotMat *  posAtom;
		// add the relative translation with is now in the coordinate system of ligand[lig]
		posAtom += tRel;

		xTr[i] = posAtom.x;
		yTr[i] = posAtom.y;
		zTr[i] = posAtom.z;
 	}
 }
>>>>>>> AS35_mb_scoring_func



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
