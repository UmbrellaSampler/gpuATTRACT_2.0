/*
 * RotMatFunctions.h
 *
 *  Created on: Aug 13, 2016
 *      Author: uwe
 */

#ifndef SRC_MATRIXFUNCTIONS_H_
#define SRC_MATRIXFUNCTIONS_H_

#include "RotMat.h"
#include "TorqueMat.h"
#include "nativeTypesWrapper.h"

#ifndef __CUDACC__

#include <cmath>
#include <ostream>
using std::cos;
using std::sin;
using std::fmin;
using std::fmax;
using std::atan2;
using std::acos;

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

#endif

namespace as {

template<typename REAL>
__host__ __device__
RotMat<REAL>  euler2rotmat(REAL const& phi, REAL const& ssi, REAL const& rot) {
	//first rotate using rot
	//		this is a rotation around the Z axis
	//			rotating Y into X and X into -Y
	//then rotate using ssi
	//		this is a rotation around the (new) Y axis
	//			rotating X into -Z and Z into X
	//finally, rotate using phi
	//		this is a rotation around the (new) Z axis
	//			rotating X into Y and Y into -X

	REAL cSSI = cos(ssi);
	REAL cPHI = cos(phi);
	REAL sSSI = sin(ssi);
	REAL sPHI = sin(phi);

	REAL cSSI_cPHI = cSSI * cPHI;
	REAL cSSI_sPHI = cSSI * sPHI;
	REAL sSSI_cPHI = sSSI * cPHI;
	REAL sSSI_sPHI = sSSI * sPHI;
	REAL cROT = cos(rot);
	REAL sROT = sin(rot);

	RotMat<REAL> rotmat;
	rotmat.mat[0] = cROT * cSSI_cPHI + sROT * sPHI;
	rotmat.mat[1] = sROT * cSSI_cPHI - cROT * sPHI;
	rotmat.mat[2] = sSSI_cPHI;

	rotmat.mat[3] = cROT * cSSI_sPHI - sROT * cPHI;
	rotmat.mat[4] = sROT * cSSI_sPHI + cROT * cPHI;
	rotmat.mat[5] = sSSI_sPHI;

	rotmat.mat[6] = -cROT * sSSI;
	rotmat.mat[7] = -sROT * sSSI;
	rotmat.mat[8] = cSSI;

	return rotmat;
}

template<typename REAL>
__host__ __device__
TorqueMat<REAL> euler2torquemat(const REAL& phi, const REAL& ssi, const REAL& rot)
{
	REAL cSSI = cos(ssi);
	REAL cPHI = cos(phi);
	REAL sSSI = sin(ssi);
	REAL sPHI = sin(phi);

	REAL cSSI_cPHI = cSSI  *  cPHI;
	REAL cSSI_sPHI = cSSI  *  sPHI;
	REAL sSSI_cPHI = sSSI  *  cPHI;
	REAL sSSI_sPHI = sSSI  *  sPHI;
	REAL cROT = cos(rot);
	REAL sROT = sin(rot);

	TorqueMat<REAL> torqueMat;
	torqueMat.mat[0][0][0] = -cROT * cSSI_sPHI + sROT * cPHI;
	torqueMat.mat[0][0][1] = -sROT * cSSI_sPHI - cROT * cPHI; // -rotmat[4]
	torqueMat.mat[0][0][2] = -sSSI_sPHI; // -rotmat[5]
	torqueMat.mat[1][0][0] = cROT * cSSI_cPHI + sROT * sPHI; // rotmat[0]
	torqueMat.mat[1][0][1] = sROT * cSSI_cPHI - cROT * sPHI; // rotmat[1]
	torqueMat.mat[1][0][2] = sSSI_cPHI; // rotmat[2]
	torqueMat.mat[2][0][0] = 0.0;
	torqueMat.mat[2][0][1] = 0.0;
	torqueMat.mat[2][0][2] = 0.0;

	torqueMat.mat[0][1][0] = -cROT * sSSI_cPHI;
	torqueMat.mat[0][1][1] = -sROT * sSSI_cPHI;
	torqueMat.mat[0][1][2] = cSSI_cPHI;
	torqueMat.mat[1][1][0] = -cROT * sSSI_sPHI;
	torqueMat.mat[1][1][1] = -sROT * sSSI_sPHI;
	torqueMat.mat[1][1][2] = cSSI_sPHI;
	torqueMat.mat[2][1][0] = -cROT * cSSI;
	torqueMat.mat[2][1][1] = -sROT * cSSI; // rotmat[7]
	torqueMat.mat[2][1][2] = -sSSI;

	torqueMat.mat[0][2][0] = -sROT * cSSI_cPHI + cROT * sPHI;
	torqueMat.mat[0][2][1] = cROT * cSSI_cPHI + sROT * sPHI;
	torqueMat.mat[0][2][2] = 0.0;
	torqueMat.mat[1][2][0] = -sROT * cSSI_sPHI - cROT * cPHI;
	torqueMat.mat[1][2][1] = cROT * cSSI_sPHI - sROT * cPHI; // rotmat[3]
	torqueMat.mat[1][2][2] = 0.0;
	torqueMat.mat[2][2][0] = sROT * sSSI;
	torqueMat.mat[2][2][1] = -cROT * sSSI; // rotmat[6]
	torqueMat.mat[2][2][2] = 0.0;

	return torqueMat;
}

template<typename T>
__host__ __device__
void rotmat2euler(const RotMat<T>& rotmat, T& phi, T& ssi, T& rot) {
	phi = atan2(rotmat[5], rotmat[2]);
	ssi = acos(fmin(fmax(rotmat[8],-1.0), 1.0));
	rot = atan2(-rotmat[7], -rotmat[6]);
	/* handel gimbal lock */
	if (abs(rotmat[8] >= 0.9999)) {
		phi = 0.0;
		if(abs(rotmat[0]) >= 0.9999) {
			if(rotmat[0] < 0.0) {
				rot = M_PI;
			} else {
				rot = 0.0;
			}

			if(rotmat[8] < 0.0) {
				ssi = M_PI;
			} else {
				ssi = 0.0;
			}

		} else {
			if(rotmat[8] < 0.0) {
				ssi = M_PI;
				rot = -acos(-rotmat[0]);
			} else {
				ssi = 0.0;
				rot = acos(rotmat[0]);
			}

		}

		if (rotmat[1] < 0) {
			rot *= -1.0;
		}
	}
}

}  // namespace as



#endif /* SRC_MATRIXFUNCTIONS_H_ */
