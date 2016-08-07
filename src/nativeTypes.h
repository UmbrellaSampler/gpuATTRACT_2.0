/*
 * nativeTypes.h
 *
 *  Created on: Apr 3, 2016
 *      Author: uwe
 */

#ifndef SRC_NATIVETYPES_H_
#define SRC_NATIVETYPES_H_

#ifdef CUDA
// from cuda toolkit
#include "vector_functions.h" // CUDA vector_types + functions (e.g. make_float4(...))

#elif defined(__GNUC__)

#include "align.h"

#ifndef __device_builtin__
#define __device_builtin__
#endif

#ifndef __cuda_builtin_vector_align8(tag, members)
#define __cuda_builtin_vector_align8(tag, members) \
struct __device_builtin__ __align__(8) tag         \
{                                                  \
    members                                        \
}
#endif

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

__cuda_builtin_vector_align8(int2, int x; int y;);
__cuda_builtin_vector_align8(uint2, unsigned int x; unsigned int y;);

struct __device_builtin__ int3
{
    int x, y, z;
};

struct __device_builtin__ uint3
{
    unsigned int x, y, z;
};

struct __device_builtin__ __builtin_align__(16) int4
{
    int x, y, z, w;
};

struct __device_builtin__ __builtin_align__(16) uint4
{
    unsigned int x, y, z, w;
};

__cuda_builtin_vector_align8(float2, float x; float y;);

struct __device_builtin__ float3
{
    float x, y, z;
};

struct __device_builtin__ __builtin_align__(16) float4
{
    float x, y, z, w;
};

struct __device_builtin__ __builtin_align__(16) double2
{
    double x, y;
};

struct __device_builtin__ double3
{
    double x, y, z;
};

struct __device_builtin__ __builtin_align__(16) double4
{
    double x, y, z, w;
};

#endif

#endif /* SRC_NATIVETYPES_H_ */
