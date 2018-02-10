/*
 * reduction_modes.h
 *
 *  Created on: Nov 22, 2017
 *      Author: glenn
 */

#ifndef REDUCTION_MODES_H_
#define REDUCTION_MODES_H_

#include "Types_6D_Modes.h"
#include "reduction.h"
#include "DeviceProtein.h"
namespace as {

template<typename REAL>
class PotForce_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
public:
	real_t E;
	vec3_t pos;
	real_t modesRec[MODES_MAX_NUMBER];
	real_t modesLig[MODES_MAX_NUMBER];
};


template<typename REAL>
void reduceModeForce(
		Vec3<REAL> const& ang,
		const REAL* forceX,
		const REAL* forceY,
		const REAL* forceZ,
		const REAL* modeX,
		const REAL* modeY,
		const REAL* modeZ,
		unsigned const& numAtoms,
		unsigned const& numModes,
		REAL* result
		)
{
	//TODO: think about passing protein to function with member "isreceptor"to determine rotation
	//rotate forces into ligand frame
	const RotMat<REAL> rotMatInv = euler2rotmat(ang.x, ang.y, ang.z).getInv();
	for( int i=0; i<numModes;i++){result[i]=0;}

	for (unsigned i = 0; i < numAtoms; ++i) {
		Vec3<REAL> forceAtom(forceX[i], forceY[i], forceZ[i]);
		forceAtom = rotMatInv * forceAtom;
		for(int mode=0;mode<numModes;mode++){
				result[mode] -= forceAtom.x*modeX[i*numModes+mode]
							  + forceAtom.y*modeY[i*numModes+mode]
							  + forceAtom.z*modeZ[i*numModes+mode];
		}
	}

}

/*
 * @bief: reduceModeForce is essentially a dot product of the force vector and the modevector
 * resulting in the effective force from the modes.
 * IMPORTANT: 1. the forces used have to be rotated such that
 * they are the corresponding coordinatesystem.
 * 2. after the forces have been reduced they have to corrected by correctModeforce
 *
 */
template<typename REAL>
void reduceModeForce(
		const REAL* forceX,
		const REAL* forceY,
		const REAL* forceZ,
		const REAL* modeX,
		const REAL* modeY,
		const REAL* modeZ,
		unsigned const& numAtoms,
		unsigned const& numModes,
		REAL* result
		)
{
	//TODO: think about passing protein to function with member "isreceptor"to determine rotation
	//rotate forces into ligand frame
	for( int i=0; i<numModes;i++){result[i]=0;}
	for (unsigned i = 0; i < numAtoms; ++i) {
		for(int mode = 0; mode < numModes; mode++){
				result[mode] -= forceX[i]*modeX[i*numModes+mode]+
								forceY[i]*modeY[i*numModes+mode]+
								forceZ[i]*modeZ[i*numModes+mode];
		}
	}
}
/**
 * @brief: this function corrects for strong mode interaction.
 * In case that high forces are acting the Proteins the mode tend
 * to intoduce too much deformation.
 * Adding an exponential term which is of higher order (harmonic not enough)
 * which is multiplied by the forceconstant for each mode corrects for this.
 * The old version uses either exp=3 or exp=4.
 * The forceconstant correspoonds to the magnitude of the modevector for each mode.
 */
template<typename REAL>
void correctModeForce(
		const REAL* modeForceConstant,
		unsigned const& numModes,
		REAL* delta
		)
{
	constexpr REAL factor = 4.0;
	constexpr int exp = 4;
	REAL counterForce;

	for(int mode = 0; mode < numModes; mode++){
		counterForce=factor*modeForceConstant[mode]*pow(delta[mode],exp);
		delta[mode]=delta[mode]+counterForce;
	}
}


#ifdef CUDA
template<typename REAL,int PROTEINTYPE, bool MODES>
void deviceReduce(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		d_Protein<REAL>* protein,
		DOF_6D_Modes<REAL>* dofs,
		REAL* xPos, REAL* yPos, REAL* zPos,
		REAL *d_fx, REAL *d_fy, REAL *d_fz,
		REAL *d_E,
		REAL *d_out,
		const cudaStream_t& stream);


/* remaining reduce part that is performed on the host */
template<typename REAL, int PROTEINTYPE, bool MODES>
void h_finalReduce(
			const unsigned& numDOFs,
			d_Protein<REAL>* protein,
			DOF_6D_Modes<REAL>* dofs,
			const REAL* deviceOut,
			Result_6D_Modes<REAL>* enGrads)
{
	unsigned dofSize = 13 ;
	if(MODES){
		dofSize += protein->numModes;
	}


	for (unsigned i = 0; i < numDOFs; ++i)
	{
		auto &enGrad = enGrads[i];
		if(PROTEINTYPE == 1){
			enGrad._6D.pos.x = deviceOut[i*dofSize + 0];
			enGrad._6D.pos.y = deviceOut[i*dofSize + 1];
			enGrad._6D.pos.z = deviceOut[i*dofSize + 2];

			for(unsigned j = 0; j < 3; ++j) {
				REAL magn2 = enGrad._6D.pos.x*enGrad._6D.pos.x
						+ enGrad._6D.pos.y*enGrad._6D.pos.y
						+ enGrad._6D.pos.z*enGrad._6D.pos.z;

				if(magn2 > static_cast<REAL>(ForceLim)) {
					enGrad._6D.pos.x *= 0.01;
					enGrad._6D.pos.y *= 0.01;
					enGrad._6D.pos.z *= 0.01;
				}
			}
			enGrad._6D.E = deviceOut[i*dofSize + 3];

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
			const auto &dof = dofs[i];
			const TorqueMat<REAL> torqueMat = euler2torquemat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z);
			Vec3<REAL> result = torqueMat.rotateReduce(torque);

			enGrad._6D.ang.x = result.x;
			enGrad._6D.ang.y = result.y;
			enGrad._6D.ang.z = result.z;
		}else{
					enGrad._6D.pos.x -= deviceOut[i*dofSize + 0];
					enGrad._6D.pos.y -= deviceOut[i*dofSize + 1];
					enGrad._6D.pos.z -= deviceOut[i*dofSize + 2];
				}

		if(MODES){
			if(PROTEINTYPE == 1){
				for(int mode=0; mode < protein->numModes; mode++){
					enGrad.modesLig[mode]=deviceOut[i*dofSize + 13 + mode];
				}
				correctModeForce(
					protein->modeForce,
					protein->numModes,
					enGrad.modesLig
					);
			}
			else{
				for(int mode=0; mode < protein->numModes; mode++){
					enGrad.modesRec[mode]=deviceOut[i*dofSize + 13 + mode];
				}

				correctModeForce(
					protein->modeForce,
					protein->numModes,
					enGrad.modesRec
					);
			}

		}

	}
}
#endif // CUDA





}//end namespace as
#endif /* REDUCTION_MODES_H_ */
