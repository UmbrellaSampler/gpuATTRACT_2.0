/*
 * transform_modes.h
 *
 *  Created on: Dec 5, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TRANSFORM_MODES_H_
#define SRC_MODEL_TRANSFORM_MODES_H_


#include "transform.h"
#include "Types_6D_Modes.h"
#include "DeviceProtein.h"

 namespace as{



 /**
  * @brief: to change from the receptor system to the ligand system
  * the translational and rotational DOF have to be "inverted".
  * This function returns an inverted DOF such that it points to the receptor in the ligandsystem
  *
  */
template<typename REAL>
const DOF_6D_Modes<REAL> invertDOF( DOF_6D_Modes<REAL> ligandDOF)
{
	DOF_6D_Modes<REAL> invertedDOF;
	Vec3<REAL> ang(0.0);
	invertedDOF._6D.ang=ligandDOF._6D.ang.inv();
	invertedDOF._6D.pos=ligandDOF._6D.pos.inv();
	const RotMat<REAL> rotMat=euler2rotmat(ligandDOF._6D.ang.x, ligandDOF._6D.ang.y, ligandDOF._6D.ang.z).getInv();
	invertedDOF._6D.pos=rotMat*invertedDOF._6D.pos;
	std::copy( ligandDOF.modesRec, ligandDOF.modesRec+MODES_MAX_NUMBER,
				invertedDOF.modesRec);
	std::copy(ligandDOF.modesLig, ligandDOF.modesLig+MODES_MAX_NUMBER,
			invertedDOF.modesLig);

return invertedDOF;
}



/**
 * @brief: In addition to the normal transform function rotate_translate_deform
 * also applies mode deformation to the coodinates.
 * Note that the deformed coordinates are also returned.
 * The deformed but not translated coordinates are important to evaluate the
 * NLForces between the receptor and the ligand
 *
 */
template<typename REAL>
void rotate_translate_deform(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		Vec3<REAL> const& displacement,
		Vec3<REAL> const& ang,
		unsigned const& numAtoms,
		unsigned const& numModes,
		REAL const* dlig,
		REAL const* xModes,
		REAL const* yModes,
		REAL const* zModes,
		REAL* xDeformed,
		REAL* yDeformed,
		REAL* zDeformed,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)
{
	const RotMat<REAL> rotMat = euler2rotmat(ang.x, ang.y, ang.z);
	 Vec3<REAL> center(0.0f) ;
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

		posAtom = rotMat*posAtom;
		posAtom += displacement;

		xTr[i] = posAtom.x;
		yTr[i] = posAtom.y;
		zTr[i] = posAtom.z;
	}
}

/**
 * Overloaded Version of the above function.
 * USECASE: for the ligand its not needed to use deformed but not translated coordinates
 */
template<typename REAL>
void rotate_translate_deform(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		Vec3<REAL> const& displacement,
		Vec3<REAL> const& ang,
		unsigned const& numAtoms,
		unsigned const& numModes,
		REAL const* dlig,
		REAL const* xModes,
		REAL const* yModes,
		REAL const* zModes,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)
{
	const RotMat<REAL> rotMat = euler2rotmat(ang.x, ang.y, ang.z);
	for (unsigned i = 0; i < numAtoms; ++i) {
		Vec3<REAL> posAtom(x[i], y[i], z[i]);
		for(int mode=0;mode<numModes;mode++){
			posAtom.x+=dlig[mode]*xModes[i*numModes+mode];
			posAtom.y+=dlig[mode]*yModes[i*numModes+mode];
			posAtom.z+=dlig[mode]*zModes[i*numModes+mode];
		}
		posAtom = rotMat*posAtom;
		posAtom += displacement;

		xTr[i] = posAtom.x;
		yTr[i] = posAtom.y;
		zTr[i] = posAtom.z;
	}
}

/*
 *
 * only applies mode deformation to an array of coordinates. Not used anymore.
 */


/*
 * This function rotates an incomming array of type REAl.
 * USECASE: since the forces action on the receptor are evaluated in the ligandframe,
 * they have to be rotated back into the global/receptor system which is what this function does
 */

template<typename REAL>
void rotate_forces(
		Vec3<REAL> const& ang,
		unsigned const& numAtoms,
		REAL* forceX,
		REAL* forceY,
		REAL* forceZ
		)
{
	const RotMat<REAL> rotMat = euler2rotmat(ang.x, ang.y, ang.z);
	 Vec3<REAL> center(0.0f) ;
	for (unsigned i = 0; i < numAtoms; ++i) {
		Vec3<REAL> forceAtom(forceX[i], forceY[i], forceZ[i]);

		forceAtom = rotMat*forceAtom;

		forceX[i] = forceAtom.x;
		forceY[i] = forceAtom.y;
		forceZ[i] = forceAtom.z;

	}
}



#ifdef CUDA

template<typename REAL, int PROTEINTYPE, bool MODES>
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
		);

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
);

#endif


}

#endif /* TRANSFORM_MODES_H_ */
