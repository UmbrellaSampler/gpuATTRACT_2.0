/*
 * reduction_modes.h
 *
 *  Created on: Nov 22, 2017
 *      Author: glenn
 */

#ifndef REDUCTION_MODES_H_
#define REDUCTION_MODES_H_

#include "Types_6D_Config.h"
#include "reduction.h"

namespace as {

template<typename REAL>
class PotForce_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
public:
	real_t E;
	vec3_t pos;
	real_t modes[MODES_LIGAND_MAX_NUMBER];
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
		forceAtom = rotMatInv*forceAtom;
		for(int mode=0;mode<numModes;mode++){
				result[mode] -= forceAtom.x*modeX[i*numModes+mode]+forceAtom.y*modeY[i*numModes+mode]+forceAtom.z*modeZ[i*numModes+mode];
		}
	}
}

/*
 * @bief: reduceModeForce is essentially a dot product of the force vector and the modevector resulting in the effective force from the modes.
 * IMPORTANT: 1. the forces used have to be rotated such that they are the corresponding coordinatesystem.
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
		for(int mode=0;mode<numModes;mode++){
				result[mode] -= forceX[i]*modeX[i*numModes+mode]+forceY[i]*modeY[i*numModes+mode]+forceZ[i]*modeZ[i*numModes+mode];
		}
	}
}
/**
 * @brief: this function corrects for strong mode interaction.
 * In case that high forces are acting the Proteins the mode tend to intoduce too much deformation.
 * Adding an exponential term which is of higher order (harmonic not enough) which is multiplied by the forceconstant for each mode corrects for this.
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
	const REAL factor=4.0;
	const int exp=4;
	REAL counterForce;

	for(int mode=0; mode<numModes; mode++){
		counterForce=factor*modeForceConstant[mode]*pow(delta[mode],exp);
		delta[mode]=delta[mode]+counterForce;
	}
}


}//end namespace as
#endif /* REDUCTION_MODES_H_ */
