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
#include "Protein.h"
#include <type_traits>


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
	invertedDOF._6D.ang=ligandDOF._6D.ang;
	invertedDOF._6D.pos=ligandDOF._6D.pos.inv() ;
	const RotMat<REAL> rotMat=euler2rotmat(ligandDOF._6D.ang.x, ligandDOF._6D.ang.y, ligandDOF._6D.ang.z).getInv();
	invertedDOF._6D.pos=rotMat*invertedDOF._6D.pos  ;
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

template<typename REAL, typename DOF_T>
__inline__ void h_deform( DOF_T const& dof, Protein<REAL> const* protein, unsigned idx_protein, unsigned idx_atom, Vec3<REAL> & posAtom,REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0)
{
	unsigned numModes = protein->numModes();
	REAL const * dlig;
	if ( idx_protein == 0){
		dlig = dof.modesRec;
	}
	else{
		dlig = dof.modesLig;
	}
	for(int mode=0;mode<numModes;mode++){
		posAtom.x += dlig[mode] * protein->xModes()[idx_atom*numModes+mode];
		posAtom.y += dlig[mode] * protein->yModes()[idx_atom*numModes+mode];
		posAtom.z += dlig[mode] * protein->zModes()[idx_atom*numModes+mode];
	}

	buffer_defoX[idx_atom] = posAtom.x;
	buffer_defoY[idx_atom] = posAtom.y;
	buffer_defoZ[idx_atom] = posAtom.z;

}

template<typename REAL, typename DOF_T>
__inline__ void h_deform( DOF_T const &dof, Protein<REAL> const* protein, unsigned idx_protein, unsigned idx_atom, Vec3<REAL> & posAtom,REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0)
{}

template<typename REAL, typename DOF_T>
__inline__ void h_rotate_translate( DOF_T const&dof,  Vec3<REAL> & posAtom, unsigned const type_protein,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0)
{
	RotMat<REAL> rotMat = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
	Vec3<REAL> translation = dof._6D.pos;
	if ( type_protein == 0)
	{
		rotMat = rotMat.getInv();
		translation = rotMat * translation.inv();
	}

	posAtom = rotMat*posAtom;
	posAtom += translation;
	}

template<typename REAL, typename DOF_T>
__inline__ void h_rotate_translate( DOF_T &dof,  Vec3<REAL> & posAtom, unsigned const type_protein,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0)
{
	RotMat<REAL> rotMat = euler2rotmat(dof.ang.x, dof.ang.y, dof.ang.z);
	if( type_protein == 0){
		rotMat =  rotMat.getInv();
	}

	posAtom = rotMat*posAtom;
	posAtom += dof.pos;
	}



template< typename REAL, typename DOF_T>
void h_DOFPos(
		Protein<REAL> const* protein,
		DOF_T const& dof,
		unsigned const type_protein,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ
		){
	REAL const* x = protein->xPos();
	REAL const* y = protein->yPos();
	REAL const* z = protein->zPos();
	for ( unsigned idx_atom = 0; idx_atom < protein->numAtoms(); ++idx_atom ){
		Vec3<REAL> posAtom(x[idx_atom], y[idx_atom], z[idx_atom]);
		h_deform( dof, protein, type_protein,  idx_atom, posAtom, buffer_defoX,
				buffer_defoY,
				buffer_defoZ);
		h_rotate_translate( dof,  posAtom, type_protein);

		buffer_trafoX[idx_atom] = posAtom.x;
		buffer_trafoY[idx_atom] = posAtom.y;
		buffer_trafoZ[idx_atom] = posAtom.z;

	}

}



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

template<typename REAL, typename DOF_T >
__device__ __forceinline__ void d_DOFPos_device(
		d_Protein<REAL> const&  protein,
		DOF_T const & dof,
		unsigned const idx,
		unsigned const type_protein,
		 REAL& buffer_defoX, REAL& buffer_defoY, REAL& buffer_defoZ,
		 REAL& buffer_trafoX, REAL& buffer_trafoY, REAL& buffer_trafoZ
		);

template<typename REAL, typename DOF_T>
 void d_DOFPos(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL> const&  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		REAL* buffer_defoX, REAL* buffer_defoY, REAL* buffer_defoZ,
		 REAL* buffer_trafoX, REAL* buffer_trafoY, REAL* buffer_trafoZ
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
