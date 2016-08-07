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

#ifndef VEC3_H_
#define VEC3_H_

#include "align.h"

#ifndef __CUDACC__

#include <cmath>
#include <ostream>
using std::max;
using std::abs;
using std::sqrt;

#endif

namespace as{


#ifndef __CUDACC__


template<typename type> class Vec3;
/** Global operators */

template<typename type>
std::ostream& operator<<(std::ostream& stream, const Vec3<type>& v) {

	stream << "[" << v.x << "; " << v.y << "; " << v.z << "]";
	return stream;
}

#endif

template<typename type>
Vec3<type> operator*(double scalar, const Vec3<type>& v) {
	return v * scalar;
}

template<typename type>
class Vec3 {
public:

	Vec3() {
	}

	Vec3(type arg) {
		x = arg;
		y = arg;
		z = arg;
	}

	Vec3(type arg0, type arg1, type arg2) {
		x = arg0;
		y = arg1;
		z = arg2;
	}

	Vec3(type args[3]) {
		x = args[0];
		y = args[1];
		z = args[2];
	}

	Vec3(const Vec3& other) {
		x = other.x;
		y = other.y;
		z = other.z;
	}

	Vec3 operator+(const Vec3& rhs) const {
		Vec3 result;
		result.x = x + rhs.x;
		result.y = y + rhs.y;
		result.z = z + rhs.z;
		return Vec3(result);
	}

	void operator+=(const Vec3& rhs) {
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
	}

	Vec3 operator-(const Vec3& rhs) const {
		Vec3 result;
		result.x = x - rhs.x;
		result.y = y - rhs.y;
		result.z = z - rhs.z;
		return Vec3(result);
	}

	void operator-=(const Vec3& rhs) {
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
	}

	template<typename REAL>
	Vec3 operator*(REAL scalar) const {
		Vec3 result;
		result.x = x * static_cast<type>(scalar);
		result.y = y * static_cast<type>(scalar);
		result.z = z * static_cast<type>(scalar);
		return result;
	}

	template<typename REAL>
	void operator*=(REAL scalar) {
		x *= static_cast<type>(scalar);
		y *= static_cast<type>(scalar);
		z *= static_cast<type>(scalar);
	}

	template<typename REAL>
	Vec3 operator/(REAL scalar) const {
		Vec3 result;
		result.x = x / static_cast<type>(scalar);
		result.y = y / static_cast<type>(scalar);
		result.z = z / static_cast<type>(scalar);
		return Vec3(result);
	}

	template<typename REAL>
	void operator/=(REAL scalar) {
		x /= static_cast<type>(scalar);
		y /= static_cast<type>(scalar);
		z /= static_cast<type>(scalar);
	}

	type MaxNorm() const {
		type norm;
		norm = max(abs(x), abs(y));
		norm = max(norm, abs(z));
		return norm;
	}

	type L2NormSquare() const {
		return x * x + y * y + z * z;
	}

	double L2Norm() const {
		return sqrt(L2NormSquare());
	}

	Vec3& operator=(const Vec3& rhs) {
		if (this != &rhs) {
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
		}
		return *this;
	}

	Vec3& operator=(type rhs) {
		x = rhs;
		y = rhs;
		z = rhs;
		return *this;
	}

	bool operator==(const Vec3& rhs) const {
		if (x != rhs.x)
			return false;
		if (y != rhs.y)
			return false;
		if (z != rhs.z)
			return false;
		return true;
	}

	bool operator!=(const Vec3& rhs) const {
		return !(*this == rhs);
	}

	~Vec3() {
	}

	type x, y, z;
};


} /* namespace asUtils */

#endif /* VEC3_H_ */
