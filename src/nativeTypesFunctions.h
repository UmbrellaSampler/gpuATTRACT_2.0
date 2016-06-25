/*
 * nativeTypesFunctions.h
 *
 *  Created on: Apr 3, 2016
 *      Author: uwe
 */

#ifndef SRC_NATIVETYPESFUNCTIONS_H_
#define SRC_NATIVETYPESFUNCTIONS_H_

#if defined(__CUDACC__)

#include "vector_functions.h"

#elif defined(__GNUC__)

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

#endif /* SRC_NATIVETYPESFUNCTIONS_H_ */
