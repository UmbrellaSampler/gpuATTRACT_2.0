/*
 * nativeTypesFunctions.h
 *
 *  Created on: Apr 3, 2016
 *      Author: uwe
 */

#ifndef SRC_NATIVETYPESFUNCTIONS_H_
#define SRC_NATIVETYPESFUNCTIONS_H_

#include "nativeTypesWrapper.h"

#ifdef CUDA

#include "vector_functions.h"

#else

#include "nativeTypes.h"

#ifndef __VECTOR_FUNCTIONS_DECL__
#define __VECTOR_FUNCTIONS_DECL__ static inline __host__ __device__
#endif

__VECTOR_FUNCTIONS_DECL__ int2 make_int2(int x, int y)
{
  int2 t; t.x = x; t.y = y; return t;
}

__VECTOR_FUNCTIONS_DECL__ uint2 make_uint2(unsigned int x, unsigned int y)
{
  uint2 t; t.x = x; t.y = y; return t;
}

__VECTOR_FUNCTIONS_DECL__ int3 make_int3(int x, int y, int z)
{
  int3 t; t.x = x; t.y = y; t.z = z; return t;
}

__VECTOR_FUNCTIONS_DECL__ uint3 make_uint3(unsigned int x, unsigned int y, unsigned int z)
{
  uint3 t; t.x = x; t.y = y; t.z = z; return t;
}

__VECTOR_FUNCTIONS_DECL__ int4 make_int4(int x, int y, int z, int w)
{
  int4 t; t.x = x; t.y = y; t.z = z; t.w = w; return t;
}

__VECTOR_FUNCTIONS_DECL__ uint4 make_uint4(unsigned int x, unsigned int y, unsigned int z, unsigned int w)
{
  uint4 t; t.x = x; t.y = y; t.z = z; t.w = w; return t;
}

__VECTOR_FUNCTIONS_DECL__ float2 make_float2(float x, float y)
{
  float2 t; t.x = x; t.y = y; return t;
}

__VECTOR_FUNCTIONS_DECL__ float3 make_float3(float x, float y, float z)
{
  float3 t; t.x = x; t.y = y; t.z = z; return t;
}

__VECTOR_FUNCTIONS_DECL__ float4 make_float4(float x, float y, float z, float w)
{
  float4 t; t.x = x; t.y = y; t.z = z; t.w = w; return t;
}

__VECTOR_FUNCTIONS_DECL__ double2 make_double2(double x, double y)
{
  double2 t; t.x = x; t.y = y; return t;
}

__VECTOR_FUNCTIONS_DECL__ double3 make_double3(double x, double y, double z)
{
  double3 t; t.x = x; t.y = y; t.z = z; return t;
}

__VECTOR_FUNCTIONS_DECL__ double4 make_double4(double x, double y, double z, double w)
{
  double4 t; t.x = x; t.y = y; t.z = z; t.w = w; return t;
}

#endif

//template<typename REAL>
//typename TypeWrapper<REAL>::real2_t make_real2(REAL x, REAL y) {};
//
//template<typename REAL>
//typename TypeWrapper<REAL>::real3_t make_real3(REAL x, REAL y, REAL z);
//
//template<typename REAL>
//typename TypeWrapper<REAL>::real4_t make_real4(REAL x, REAL y, REAL z, REAL w);

static inline float2 make_real2(float x, float y) {
	return make_float2(x, y);
}

static inline float3 make_real3(float x, float y, float z) {
	return make_float3(x, y, z);
}

static inline float4 make_real4(float x, float y, float z, float w) {
	return make_float4(x, y, z, w);
}

static inline double2 make_real2(double x, double y) {
	return make_double2(x, y);
}

static inline double3 make_real3(double x, double y, double z) {
	return make_double3(x, y, z);
}

static inline double4 make_real4(double x, double y, double z, double w) {
	return make_double4(x, y, z, w);
}



#endif /* SRC_NATIVETYPESFUNCTIONS_H_ */
