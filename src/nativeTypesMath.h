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

#ifndef NATIVETYPESMATH_H_
#define NATIVETYPESMATH_H_

#include "nativeTypes.h"
#include "nativeTypesFunctions.h"

inline __host__  __device__ float4 operator-(float4 a, float4 b) {
	return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

inline __host__  __device__ float4 operator+(float4 a, float4 b) {
	return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline __host__  __device__ float4 operator*(float4 a, float4 b) {
	return make_float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

inline __host__  __device__ float4 operator*(float4 a, float b) {
	return make_float4(a.x * b, a.y * b, a.z * b, a.w * b);
}

inline __host__  __device__ float4 operator/(float4 a, float4 b) {
	return make_float4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

inline __host__  __device__ float4 operator/(float4 a, float b) {
	return make_float4(a.x / b, a.y / b, a.z / b, a.w / b);
}

inline __host__  __device__ float4 operator-(float4 a) {
	return make_float4(-a.x , -a.y, -a.z, -a.w);
}

inline __host__  __device__ float3 operator-(float3 a, float3 b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__  __device__ float3 operator+(float3 a, float3 b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __host__  __device__ float3 operator*(float3 a, float3 b) {
	return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __host__  __device__ float3 operator*(float3 a, float b) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __host__  __device__ float3 operator/(float3 a, float3 b) {
	return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline __host__  __device__ float3 operator/(float3 a, float b) {
	return make_float3(a.x / b, a.y / b, a.z / b);
}

inline __host__  __device__ float3 operator-(float3 a) {
	return make_float3(-a.x, -a.y, -a.z);
}

/* double */

inline __host__  __device__ double4 operator-(double4 a, double4 b) {
	return make_double4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

inline __host__  __device__ double4 operator+(double4 a, double4 b) {
	return make_double4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline __host__  __device__ double4 operator*(double4 a, double4 b) {
	return make_double4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

inline __host__  __device__ double4 operator*(double4 a, double b) {
	return make_double4(a.x * b, a.y * b, a.z * b, a.w * b);
}

inline __host__  __device__ double4 operator/(double4 a, double4 b) {
	return make_double4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

inline __host__  __device__ double4 operator/(double4 a, double b) {
	return make_double4(a.x / b, a.y / b, a.z / b, a.w / b);
}

inline __host__  __device__ double4 operator-(double4 a) {
	return make_double4(-a.x, -a.y, -a.z, -a.w);
}

inline __host__  __device__ double3 operator-(double3 a, double3 b) {
	return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__  __device__ double3 operator+(double3 a, double3 b) {
	return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __host__  __device__ double3 operator*(double3 a, double3 b) {
	return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __host__  __device__ double3 operator*(double3 a, double b) {
	return make_double3(a.x * b, a.y * b, a.z * b);
}

inline __host__  __device__ double3 operator/(double3 a, double3 b) {
	return make_double3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline __host__  __device__ double3 operator/(double3 a, double b) {
	return make_double3(a.x / b, a.y / b, a.z / b);
}

inline __host__  __device__ double3 operator-(double3 a) {
	return make_double3(-a.x, -a.y, -a.z);
}


#endif /* NATIVETYPESMATH_H_ */
