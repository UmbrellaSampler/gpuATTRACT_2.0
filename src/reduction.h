/*
 * reduction.h
 *
 *  Created on: Aug 13, 2016
 *      Author: uwe
 */

#ifndef SRC_REDUCTION_H_
#define SRC_REDUCTION_H_

#include "Vec3.h"
#include "Torque.h"
#include "TorqueMat.h"
#include "MatrixFunctions.h"
#include "nativeTypesWrapper.h"

namespace as {

template<typename REAL>
class PotForce {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
public:
	real_t E;
	vec3_t pos;
};


constexpr float ForceLim = 1.0e18;

template<typename REAL>
PotForce<REAL> reducePotForce(
		REAL const* fx, REAL const* fy, REAL const* fz,
		REAL const* energy,
		unsigned const& numAtoms)
{
	PotForce<REAL> potForce;
	for(unsigned i = 0; i < numAtoms; ++i) {
		potForce.pos.x += fx[i];
		potForce.pos.y += fy[i];
		potForce.pos.z += fz[i];
		potForce.E += energy[i];
	}

	// force reduction, some times helps in case of very "bad" start structure
	// taken from original ATTRACT code in trans.f
	for(unsigned i = 0; i < 3; ++i) {
		REAL magn2 = potForce.pos.x*potForce.pos.x
				+ potForce.pos.y*potForce.pos.y
				+ potForce.pos.z*potForce.pos.z;

		if(magn2 > static_cast<REAL>(ForceLim)) {
			potForce.pos.x *= 0.01;
			potForce.pos.y *= 0.01;
			potForce.pos.z *= 0.01;
		}
	}

	return potForce;
}

template<typename REAL>
Vec3<REAL> reduceTorque(
		REAL const* x, REAL const* y, REAL const* z, // unrotated protein coordinates x,y,z !!!
		REAL const* fx, REAL const* fy, REAL const* fz,
		unsigned const& numAtoms,
		Vec3<REAL> const& ang)
{
	Torque<REAL> torque(0);
	for (unsigned i = 0; i < numAtoms; ++i) {
		torque.mat[0][0] += x[i]*fx[i];
		torque.mat[0][1] += y[i]*fx[i];
		torque.mat[0][2] += z[i]*fx[i];
		torque.mat[1][0] += x[i]*fy[i];
		torque.mat[1][1] += y[i]*fy[i];
		torque.mat[1][2] += z[i]*fy[i];
		torque.mat[2][0] += x[i]*fz[i];
		torque.mat[2][1] += y[i]*fz[i];
		torque.mat[2][2] += z[i]*fz[i];
	}

	TorqueMat<REAL> torqueMat = euler2torquemat(ang.x, ang.y, ang.z);
	Vec3<REAL> result = torqueMat.rotateReduce(torque);

	return result;
}

} // namespace as


#endif /* SRC_REDUCTION_H_ */
