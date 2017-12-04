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
		forceAtom = rotMatInv*forceAtom;
		for(int mode=0;mode<numModes;mode++){
				result[mode] -= forceAtom.x*modeX[i*numModes+mode]+forceAtom.y*modeY[i*numModes+mode]+forceAtom.z*modeZ[i*numModes+mode];
		}
	}

}

}//end namespace as
#endif /* REDUCTION_MODES_H_ */
