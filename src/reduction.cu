/*
 * reduction.cu
 *
 *  Created on: Sep 5, 2016
 *      Author: uwe
 */

#include "reduction.h"

namespace as {

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator       double *()
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};

// each thread block reduces the arrays one structure
template <typename T, unsigned int blockSize, bool nIsPow2>
__global__ void
blockReduce(unsigned numAtoms,
		T* xPos, T* yPos, T* zPos,
		T *d_fx, T *d_fy, T *d_fz, T *d_E,
		T *g_odata)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = threadIdx.x; // half the number of blocks
    unsigned int base = blockIdx.x*numAtoms;

    T sum_fx = 0;
    T sum_fy = 0;
    T sum_fz = 0;
    T sum_E = 0;
    T sum_torque[9] = {0};

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < numAtoms)
    {
    	T fx, fy, fz, x, y, z;
		fx = d_fx[base + i];
		fy = d_fy[base + i];
		fz = d_fz[base + i];
		x = xPos[i];
		y = yPos[i];
		z = zPos[i];
		sum_fx += fx;
		sum_fy += fy;
		sum_fz += fz;
		sum_E += d_E[base + i];
		sum_torque[0] += x * fx;
		sum_torque[1] += y * fx;
		sum_torque[2] += z * fx;
		sum_torque[3] += x * fy;
		sum_torque[4] += y * fy;
		sum_torque[5] += z * fy;
		sum_torque[6] += x * fz;
		sum_torque[7] += y * fz;
		sum_torque[8] += z * fz;


        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < numAtoms) {

			fx = d_fx[base + i + blockSize];
			fy = d_fy[base + i + blockSize];
			fz = d_fz[base + i + blockSize];
			x = xPos[i + blockSize];
			y = yPos[i + blockSize];
			z = zPos[i + blockSize];
			sum_fx += fx;
			sum_fy += fy;
			sum_fz += fz;
			sum_E += d_E[base + i + blockSize];
			sum_torque[0] += x * fx;
			sum_torque[1] += y * fx;
			sum_torque[2] += z * fx;
			sum_torque[3] += x * fy;
			sum_torque[4] += y * fy;
			sum_torque[5] += z * fy;
			sum_torque[6] += x * fz;
			sum_torque[7] += y * fz;
			sum_torque[8] += z * fz;
		}

        i += blockSize*2;
    }

    // each thread puts its local sum into shared memory
    sdata[tid + 0 * blockSize] = sum_fx;
    sdata[tid + 1 * blockSize] = sum_fy;
    sdata[tid + 2 * blockSize] = sum_fz;
    sdata[tid + 3 * blockSize] = sum_E;
    sdata[tid + 4 * blockSize] = sum_torque[0];
    sdata[tid + 5 * blockSize] = sum_torque[1];
    sdata[tid + 6 * blockSize] = sum_torque[2];
    sdata[tid + 7 * blockSize] = sum_torque[3];
    sdata[tid + 8 * blockSize] = sum_torque[4];
    sdata[tid + 9 * blockSize] = sum_torque[5];
    sdata[tid + 10* blockSize] = sum_torque[6];
    sdata[tid + 11* blockSize] = sum_torque[7];
    sdata[tid + 12* blockSize] = sum_torque[8];

    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 1024) && (tid < 512))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 512];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 512];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 512];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 512];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 512];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 512];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 512];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 512];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 512];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 512];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 512];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 512];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 512];
    }

    __syncthreads();

    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 256];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 256];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 256];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 256];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 256];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 256];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 256];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 256];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 256];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 256];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 256];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 256];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 256];
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
		sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 128];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 128];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 128];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 128];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 128];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 128];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 128];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 128];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 128];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 128];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 128];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 128];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 128];
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 64];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 64];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 64];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 64];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 64];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 64];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 64];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 64];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 64];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 64];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 64];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 64];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 64];
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) {
        	sum_fx        += sdata[tid + 0 * blockSize + 32];
			sum_fy        += sdata[tid + 1 * blockSize + 32];
			sum_fz        += sdata[tid + 2 * blockSize + 32];
			sum_E         += sdata[tid + 3 * blockSize + 32];
			sum_torque[0] += sdata[tid + 4 * blockSize + 32];
			sum_torque[1] += sdata[tid + 5 * blockSize + 32];
			sum_torque[2] += sdata[tid + 6 * blockSize + 32];
			sum_torque[3] += sdata[tid + 7 * blockSize + 32];
			sum_torque[4] += sdata[tid + 8 * blockSize + 32];
			sum_torque[5] += sdata[tid + 9 * blockSize + 32];
			sum_torque[6] += sdata[tid + 10* blockSize + 32];
			sum_torque[7] += sdata[tid + 11* blockSize + 32];
			sum_torque[8] += sdata[tid + 12* blockSize + 32];

        }
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2)
        {
            sum_fx        += __shfl_down(sum_fx       , offset);
			sum_fy        += __shfl_down(sum_fy       , offset);
			sum_fz        += __shfl_down(sum_fz       , offset);
			sum_E         += __shfl_down(sum_E        , offset);
			sum_torque[0] += __shfl_down(sum_torque[0], offset);
			sum_torque[1] += __shfl_down(sum_torque[1], offset);
			sum_torque[2] += __shfl_down(sum_torque[2], offset);
			sum_torque[3] += __shfl_down(sum_torque[3], offset);
			sum_torque[4] += __shfl_down(sum_torque[4], offset);
			sum_torque[5] += __shfl_down(sum_torque[5], offset);
			sum_torque[6] += __shfl_down(sum_torque[6], offset);
			sum_torque[7] += __shfl_down(sum_torque[7], offset);
			sum_torque[8] += __shfl_down(sum_torque[8], offset);
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 32];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 32];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 32];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 32];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 32];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 32];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 32];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 32];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 32];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 32];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 32];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 32];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 32];
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 16];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 16];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 16];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 16];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 16];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 16];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 16];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 16];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 16];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 16];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 16];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 16];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 16];
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 8];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 8];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 8];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 8];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 8];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 8];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 8];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 8];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 8];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 8];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 8];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 8];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 8];
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 4];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 4];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 4];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 4];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 4];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 4];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 4];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 4];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 4];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 4];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 4];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 4];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 4];
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 2];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 2];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 2];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 2];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 2];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 2];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 2];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 2];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 2];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 2];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 2];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 2];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 2];
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid + 0 * blockSize] = sum_fx        = sum_fx        + sdata[tid + 0 * blockSize + 1];
		sdata[tid + 1 * blockSize] = sum_fy        = sum_fy        + sdata[tid + 1 * blockSize + 1];
		sdata[tid + 2 * blockSize] = sum_fz        = sum_fz        + sdata[tid + 2 * blockSize + 1];
		sdata[tid + 3 * blockSize] = sum_E         = sum_E         + sdata[tid + 3 * blockSize + 1];
		sdata[tid + 4 * blockSize] = sum_torque[0] = sum_torque[0] + sdata[tid + 4 * blockSize + 1];
		sdata[tid + 5 * blockSize] = sum_torque[1] = sum_torque[1] + sdata[tid + 5 * blockSize + 1];
		sdata[tid + 6 * blockSize] = sum_torque[2] = sum_torque[2] + sdata[tid + 6 * blockSize + 1];
		sdata[tid + 7 * blockSize] = sum_torque[3] = sum_torque[3] + sdata[tid + 7 * blockSize + 1];
		sdata[tid + 8 * blockSize] = sum_torque[4] = sum_torque[4] + sdata[tid + 8 * blockSize + 1];
		sdata[tid + 9 * blockSize] = sum_torque[5] = sum_torque[5] + sdata[tid + 9 * blockSize + 1];
		sdata[tid + 10* blockSize] = sum_torque[6] = sum_torque[6] + sdata[tid + 10* blockSize + 1];
		sdata[tid + 11* blockSize] = sum_torque[7] = sum_torque[7] + sdata[tid + 11* blockSize + 1];
		sdata[tid + 12* blockSize] = sum_torque[8] = sum_torque[8] + sdata[tid + 12* blockSize + 1];
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) {
    	g_odata[0  + blockIdx.x*14] = sum_fx;
    	g_odata[1  + blockIdx.x*14] = sum_fy;
    	g_odata[2  + blockIdx.x*14] = sum_fz;
    	g_odata[3  + blockIdx.x*14] = sum_E;
    	g_odata[4  + blockIdx.x*14] = sum_torque[0];
    	g_odata[5  + blockIdx.x*14] = sum_torque[1];
    	g_odata[6  + blockIdx.x*14] = sum_torque[2];
    	g_odata[7  + blockIdx.x*14] = sum_torque[3];
    	g_odata[8  + blockIdx.x*14] = sum_torque[4];
    	g_odata[9  + blockIdx.x*14] = sum_torque[5];
    	g_odata[10 + blockIdx.x*14] = sum_torque[6];
    	g_odata[11 + blockIdx.x*14] = sum_torque[7];
    	g_odata[12 + blockIdx.x*14] = sum_torque[8];
    }
}

template <class T>
void d_reduce(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtoms,
		T* xPos, T* yPos, T* zPos,
		T *d_fx, T *d_fy, T *d_fz, T *d_E,
		T *g_odata,
		const cudaStream_t& stream)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize = (threads <= 32) ? 2 * 13 * threads * sizeof(T) : 13 * threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch

	if (isPow2(numAtoms))
	{
		switch (threads)
		{
			case 1024:
				blockReduce<T, 1024,true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 512:
				blockReduce<T, 512, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 256:
				blockReduce<T, 256, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 128:
				blockReduce<T, 128, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 64:
				blockReduce<T,  64, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 32:
				blockReduce<T,  32, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 16:
				blockReduce<T,  16, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  8:
				blockReduce<T,   8, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  4:
				blockReduce<T,   4, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  2:
				blockReduce<T,   2, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  1:
				blockReduce<T,   1, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;
		}
	}
	else
	{
		switch (threads)
		{
			case 1024:
				blockReduce<T, 1024,false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 512:
				blockReduce<T, 512, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 256:
				blockReduce<T, 256, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 128:
				blockReduce<T, 128, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 64:
				blockReduce<T,  64, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 32:
				blockReduce<T,  32, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 16:
				blockReduce<T,  16, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  8:
				blockReduce<T,   8, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  4:
				blockReduce<T,   4, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  2:
				blockReduce<T,   2, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  1:
				blockReduce<T,   1, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, xPos, yPos, zPos, d_fx, d_fy, d_fz, d_E, g_odata);
				break;
		}
	}

}

template
void d_reduce<float>(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtoms,
		float* xPos, float* yPos, float* zPos,
		float *d_fx, float *d_fy, float *d_fz, float *d_E,
		float *g_odata,
		const cudaStream_t& stream);

template
void d_reduce<double>(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtoms,
		double* xPos, double* yPos, double* zPos,
		double *d_fx, double *d_fy, double *d_fz, double *d_E,
		double *g_odata,
		const cudaStream_t& stream);

}  // namespace as


