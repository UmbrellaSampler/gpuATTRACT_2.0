/*
 * reduction.cu
 *
 *  Created on: Sep 5, 2016
 *      Author: uwe
 */

#include "reduction_modes.h"

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

///BLOCK REDUCE MODE////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// each thread block reduces the arrays one structure
template <typename T, unsigned int blockSize, bool nIsPow2>
__global__ void
blockReduceMode(unsigned numAtomsRec, unsigned numAtomsLig, unsigned numModesRec, unsigned numModesLig, typename Types_6D_modes<T>::DOF* dofs,
		T* xPos, T* yPos, T* zPos,
		T *xModesRec,T *yModesRec,T *zModesRec,
		T *xModesLig,T *yModesLig,T *zModesLig,
		T *d_fxRec, T *d_fyRec, T *d_fzRec,
		T *d_fxLig, T *d_fyLig, T *d_fzLig,
		T *d_E,
		T *g_odata)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = threadIdx.x; // half the number of blocks
    unsigned int maxAtom = max(numAtomsLig, numAtomsRec);
    unsigned int base = blockIdx.x* maxAtom;

    T sum_fx = 0;
    T sum_fy = 0;
    T sum_fz = 0;
    T sum_E = 0;
    T sum_torque[9] = {0};
    T sum_modeLig[10] = {0};
    T sum_modeRec[10] = {0};

    unsigned DOFidx = tid / maxAtom;
    auto dof = dofs[DOFidx];
    Vec3<T> const& ang = dof._6D.ang;
    const RotMat<T> rotMatInv = euler2rotmat(ang.x, ang.y, ang.z).getInv();

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < maxAtom)
    {

    	if(i < numAtomsRec){
			//reduce mode force for the receptor
			T fxRec, fyRec, fzRec;
			fxRec = d_fxRec[base + i];
			fyRec = d_fyRec[base + i];
			fzRec = d_fzRec[base + i];
			Vec3<T> forceAtomRec(fxRec, fyRec, fzRec);
			for(int mode=0;mode<numModesRec;mode++){
				sum_modeRec[mode] -=  forceAtomRec.x*xModesLig[i*numModesLig+mode]+forceAtomRec.y*yModesLig[i*numModesLig+mode]+forceAtomRec.z*zModesLig[i*numModesLig+mode];
			}

			if (nIsPow2 || i + blockSize < maxAtom) {
				for(int mode=0;mode<numModesRec;mode++){
					sum_modeRec[mode] +=  forceAtomRec.x*xModesRec[i*numModesRec+mode]+forceAtomRec.y*yModesRec[i*numModesRec+mode]+forceAtomRec.z*zModesRec[i*numModesRec+mode];
				}
			}
		}


    	if(i < numAtomsLig){
			T fxLig, fyLig, fzLig, x, y, z;
			fxLig = d_fxLig[base + i];
			fyLig = d_fyLig[base + i];
			fzLig = d_fzLig[base + i];
			x = xPos[i];
			y = yPos[i];
			z = zPos[i];
			sum_fx += fxLig;
			sum_fy += fyLig;
			sum_fz += fzLig;
			sum_E += d_E[base + i];
			sum_torque[0] += x * fxLig;
			sum_torque[1] += y * fxLig;
			sum_torque[2] += z * fxLig;
			sum_torque[3] += x * fyLig;
			sum_torque[4] += y * fyLig;
			sum_torque[5] += z * fyLig;
			sum_torque[6] += x * fzLig;
			sum_torque[7] += y * fzLig;
			sum_torque[8] += z * fzLig;

			//reduce mode Force for the ligand
			Vec3<T> forceAtomLig(fxLig, fyLig, fzLig);
			forceAtomLig = rotMatInv*forceAtomLig;
			for(int mode=0;mode<numModesLig;mode++){
				sum_modeLig[mode] -=  forceAtomLig.x*xModesLig[i*numModesLig+mode]+forceAtomLig.y*yModesLig[i*numModesLig+mode]+forceAtomLig.z*zModesLig[i*numModesLig+mode];
			}



        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < maxAtom) {

				fxLig = d_fxLig[base + i + blockSize];
				fyLig = d_fyLig[base + i + blockSize];
				fzLig = d_fzLig[base + i + blockSize];
				x = xPos[i + blockSize];
				y = yPos[i + blockSize];
				z = zPos[i + blockSize];
				sum_fx += fxLig;
				sum_fy += fyLig;
				sum_fz += fzLig;
				sum_E += d_E[base + i + blockSize];
				sum_torque[0] += x * fxLig;
				sum_torque[1] += y * fxLig;
				sum_torque[2] += z * fxLig;
				sum_torque[3] += x * fy;
				sum_torque[4] += y * fyLig;
				sum_torque[5] += z * fyLig;
				sum_torque[6] += x * fzLig;
				sum_torque[7] += y * fzLig;
				sum_torque[8] += z * fzLig;
				for(int mode=0;mode<numModesLig;mode++){
					sum_modeLig[mode] +=  forceAtomLig.x*xModesLig[i*numModesLig+mode]+forceAtomLig.y*yModesLig[i*numModesLig+mode]+forceAtomLig.z*zModesLig[i*numModesLig+mode];
				}
        	}

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
    for(int mode=0;mode<numModesLig;mode++){
    	sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode];
    }
    for(int mode=0;mode<numModesRec;mode++){
		sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode];
	}
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
		for(int mode=0;mode<numModesLig;mode++){
		sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 512];
    	}
    	for(int mode=0;mode<numModesRec;mode++){
    		sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 512];
        }
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13 + mode)* blockSize + 256];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 256];
		}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 128];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 128];
		}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13 + mode)* blockSize + 64];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 64];
		}
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
			for(int mode=0;mode<numModesLig;mode++){
				sum_modeLig[mode] += sdata[tid + (13 + mode)* blockSize + 32];
			}
			for(int mode=0;mode<numModesRec;mode++){
				sum_modeRec[mode] += sdata[tid + (13 + numModesLig + mode)* blockSize + 32];
			}

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
			for(int mode=0;mode<numModesLig;mode++){
				sum_modeLig[mode] +=__shfl_down(sum_modeLig[mode], offset);
			}
			for(int mode=0;mode<numModesRec;mode++){
				sum_modeRec[mode] +=__shfl_down(sum_modeRec[mode], offset);
			}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 32];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 32];
		}

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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 16];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 16];
		}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 8];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 8];
		}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 4];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 4];
		}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 2];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 2];
		}
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
		for(int mode=0;mode<numModesLig;mode++){
			sdata[tid + (13+mode)* blockSize] =  sum_modeLig[mode]=sum_modeLig[mode] + sdata[tid + (13+mode)* blockSize + 1];
		}
		for(int mode=0;mode<numModesRec;mode++){
			sdata[tid + (13+numModesLig+mode)* blockSize] =  sum_modeRec[mode]=sum_modeRec[mode] + sdata[tid + (13 + numModesLig + mode)* blockSize + 1];
		}
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) {
    	g_odata[0  + blockIdx.x*(13+numModesLig + numModesRec + numModesRec)] = sum_fx;
    	g_odata[1  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_fy;
    	g_odata[2  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_fz;
    	g_odata[3  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_E;
    	g_odata[4  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[0];
    	g_odata[5  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[1];
    	g_odata[6  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[2];
    	g_odata[7  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[3];
    	g_odata[8  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[4];
    	g_odata[9  + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[5];
    	g_odata[10 + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[6];
    	g_odata[11 + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[7];
    	g_odata[12 + blockIdx.x*(13+numModesLig + numModesRec)] = sum_torque[8];
    	for(int mode=0;mode<numModesLig;mode++){
    		g_odata[13 + mode + blockIdx.x*(13+numModesLig + numModesRec)] =  sum_modeLig[mode];
    	}
    	for(int mode=0;mode<numModesRec;mode++){
			g_odata[13 + numModesLig + mode + blockIdx.x*(13+numModesRec+numModesLig)] =  sum_modeRec[mode];
		}
    }
}
///END BLOCK REDUCE MODE



template
void d_reduce<float>(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtomsRec,
		float* xPos, float* yPos, float* zPos,
		float *d_fx, float *d_fy, float *d_fz, float *d_E,
		float *g_odata,
		const cudaStream_t& stream);

template
void d_reduce<double>(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtomsRec,
		double* xPos, double* yPos, double* zPos,
		double *d_fx, double *d_fy, double *d_fz, double *d_E,
		double *g_odata,
		const cudaStream_t& stream);
//GLENN add



template <class T>
void d_reduceMode(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtomsRec,
		const unsigned& numAtomsLig,
		const unsigned& numModesRec,
		const unsigned& numModesLig,
		typename Types_6D_modes<T>::DOF* dofs,
		T* xPos, T* yPos, T* zPos,
		T *xModesLig,T *yModesLig,T *zModesLig,
		T *xModesRec,T *yModesRec,T *zModesRec,
		T *d_fxRec, T *d_fyRec, T *d_fzRec,
		T *d_fx, T *d_fy, T *d_fz, T *d_E,
		T *g_odata,
		const cudaStream_t& stream)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	const int smemSize = (threads <= 32) ? 2 * (13 + numModesLig + numModesRec) * threads * sizeof(T) : (13 + numModesLig + numModesRec) * threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch




	if (isPow2(numAtoms))
	{
		switch (threads)
		{
			case 1024:
				blockReduceMode<T, 1024,true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 512:
				blockReduceMode<T, 512, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 256:
				blockReduceMode<T, 256, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 128:
				blockReduceMode<T, 128, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 64:
				blockReduceMode<T,  64, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 32:
				blockReduceMode<T,  32, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 16:
				blockReduceMode<T,  16, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  8:
				blockReduceMode<T,   8, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  4:
				blockReduceMode<T,   4, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  2:
				blockReduceMode<T,   2, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  1:
				blockReduceMode<T,   1, true><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec, numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec, d_fx, d_fy, d_fz, d_E, g_odata);
				break;
		}
	}
	else
	{
		switch (threads)
		{
			case 1024:
				blockReduceMode<T, 1024,false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 512:
				blockReduceMode<T, 512, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 256:
				blockReduceMode<T, 256, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 128:
				blockReduceMode<T, 128, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 64:
				blockReduceMode<T,  64, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 32:
				blockReduceMode<T,  32, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case 16:
				blockReduceMode<T,  16, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  8:
				blockReduceMode<T,   8, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  4:
				blockReduceMode<T,   4, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  2:
				blockReduceMode<T,   2, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;

			case  1:
				blockReduceMode<T,   1, false><<< dimGrid, dimBlock, smemSize, stream >>>(numAtomsRec, numAtomsLig, numModesRec,  numModesLig, dofs, xPos, yPos, zPos, xModesRec, yModesRec, zModesRec,  xModesLig, yModesLig, zModesLig, d_fxRec, d_fyRec, d_fzRec,d_fx, d_fy, d_fz, d_E, g_odata);
				break;
		}
	}

}



template
void d_reduceMode<float>(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtomsRec,
		const unsigned& numAtomsLig,
		const unsigned& numModesRec,
		const unsigned& numModesLig,
		typename Types_6D_modes<float>::DOF* dofs,
		float* xPos, float* yPos, float* zPos,
		float *xModesLig,float *yModesLig,float *zModesLig,
		float *xModesRec,float *yModesRec,float *zModesRec,
		float *d_fxRec, float *d_fyRec, float *d_fzRec,
		float *d_fx, float *d_fy, float *d_fz, float *d_E,
		float *g_odata,
		const cudaStream_t& stream);

template
void d_reduceMode<double>(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtomsRec,
		const unsigned& numAtomsLig,
		const unsigned& numModesRec,
		const unsigned& numModesLig,
		typename Types_6D_modes<double>::DOF* dofs,
		double* xPos, double* yPos, double* zPos,
		double *xModesLig,double *yModesLig,double *zModesLig,
		double *xModesRec,double *yModesRec,double *zModesRec,
		double *d_fxRec, double *d_fyRec, double *d_fzRec,
		double *d_fx, double *d_fy, double *d_fz, double *d_E,
		double *g_odata,
		const cudaStream_t& stream);


//GLENN ADD MODEFORCE


}  // namespace as


