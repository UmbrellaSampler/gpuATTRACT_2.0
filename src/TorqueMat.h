/*
 * TorqueMat.h
 *
 *  Created on: Aug 13, 2016
 *      Author: uwe
 */

#ifndef SRC_TORQUEMAT_H_
#define SRC_TORQUEMAT_H_

#include "Torque.h"
#include "Vec3.h"

namespace as {

template<typename T>
class TorqueMat {
public:
	TorqueMat() {};
	TorqueMat(T value) {
		mat = {value};
	}

	/**
	 * Rotation and Reduction of torque components
	 */
	Vec3<T> rotateReduce(Torque<T> const& torque) {
		Vec3<T> result(0.0);
		for (unsigned k = 0; k < 3; ++k) {
			for (unsigned l = 0; l < 3; ++l) {
				result.x += mat[k][0][l] * torque.mat[k][l];
				result.y += mat[k][1][l] * torque.mat[k][l];
				result.z += mat[k][2][l] * torque.mat[k][l];
			}
		}
		return result;
	}


	T mat[3][3][3];
};

} // namespace as



#endif /* SRC_TORQUEMAT_H_ */
