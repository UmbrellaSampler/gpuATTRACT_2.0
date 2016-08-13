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

#ifndef __CUDACC__

#include <cmath>
#include <ostream>
using std::cos;
using std::sin;

#endif

namespace as {

template<typename REAL>
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





}  // namespace as



#endif /* SRC_MATRIXFUNCTIONS_H_ */
