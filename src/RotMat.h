/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef ROTMAT_H_
#define ROTMAT_H_

#include <Vec3.h>
#include <algorithm>

#ifndef __CUDACC__

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

#endif

namespace as {

#ifndef __CUDACC__
template <typename T> class RotMat;
// ostream is not available in nvcc

template<typename type>
std::ostream& operator<<(std::ostream& stream, const RotMat<type>& v) {
	stream << "[";
	unsigned j = 0;
	for(int i = 0; i < 3; ++i) {
		stream << "[" << v.mat[j+0] << ", " << v.mat[j+1] << ", " << v.mat[j+2] << "]";
		j += 3;
	}
	stream << "]";
	return stream;
}
#endif

template <typename T>
class RotMat {
public:
	__host__ __device__
	RotMat() {};

	__host__ __device__
	RotMat(T value) {
		std::fill_n(mat, 9, value);
	}

	T mat[9];

	__host__ __device__
	RotMat getInv() const {
		RotMat m0;
		m0.mat[0] = mat[0];
		m0.mat[1] = mat[3];
		m0.mat[2] = mat[6];
		m0.mat[3] = mat[1];
		m0.mat[4] = mat[4];
		m0.mat[5] = mat[7];
		m0.mat[6] = mat[2];
		m0.mat[7] = mat[5];
		m0.mat[8] = mat[8];
		return m0;
	}

	__host__ __device__
	RotMat operator* (const RotMat& rhs) {
		RotMat matOut(0.0);
		for(unsigned i = 0; i < 3; ++i) {
			for(unsigned j = 0; j < 3; ++j) {
				for(unsigned k = 0; k < 3; ++k) {
					matOut[i*3 + j] += mat[i*3 + k]*rhs[k*3 + j];
				}
			}
		}
		return matOut;
	}

	__host__ __device__
	T& operator[](const unsigned& i) {
		return mat[i];
	}

	__host__ __device__
	const T& operator[](const unsigned& i) const {
		return mat[i];
	}

	__host__ __device__
	void operator*= (const RotMat& rhs) {
		RotMat matOut(0.0);
		for(unsigned i = 0; i < 3; ++i) {
			for(unsigned j = 0; j < 3; ++j) {
				for(unsigned k = 0; k < 3; ++k) {
					matOut.mat[i*3 + j] += mat[i*3 + k]*rhs.mat[k*3 + j];
				}
			}
		}
		*this = matOut;
	}

	__host__ __device__
	Vec3<T> operator * (const Vec3<T>& rhs) const noexcept {
		Vec3<T> vecOut(0.0);
		vecOut.x = mat[0] * rhs.x + mat[1] * rhs.y + mat[2] * rhs.z;
		vecOut.y = mat[3] * rhs.x + mat[4] * rhs.y + mat[5] * rhs.z;
		vecOut.z = mat[6] * rhs.x + mat[7] * rhs.y + mat[8] * rhs.z;
		return vecOut;
	}

};


}
#endif /* ROTMAT_H_ */
