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

template <typename REAL, typename DOF_T, int PROTEIN_T>
__device__ inline void rotateForce( DOF_T const & dof, REAL* forceAtom_in, REAL* forceAtom_out, unsigned const idx_protein,
		typename std::enable_if< std::is_same<DOF_T,DOF_MB_Modes<REAL>>::value, void>::type* dummy = 0)
		//typename std::enable_if< std::is_same<DOF_T,DOF_6D_Modes<REAL>>, void>::type* dummy1 = 0)
		{
	if(PROTEIN_T == 0){
				forceAtom_out[0] = forceAtom_in[0];
				forceAtom_out[1] = forceAtom_in[1];
				forceAtom_out[2] = forceAtom_in[2];
			}
	else{
		RotMat<REAL> const rotMatInv = euler2rotmat(dof.protein[idx_protein].ang.x, dof.protein[idx_protein].ang.y, dof.protein[idx_protein].ang.z).getInv();
		forceAtom_out[0] = rotMatInv[0] * forceAtom_in[0] + rotMatInv[1] * forceAtom_in[1] + rotMatInv[2] * forceAtom_in[2];
		forceAtom_out[1] = rotMatInv[3] * forceAtom_in[0] + rotMatInv[4] * forceAtom_in[1] + rotMatInv[5] * forceAtom_in[2];
		forceAtom_out[2] = rotMatInv[6] * forceAtom_in[0] + rotMatInv[7] * forceAtom_in[1] + rotMatInv[8] * forceAtom_in[2];
	}
}

template <typename REAL, typename DOF_T, int PROTEIN_T>
__device__ inline void rotateForce( DOF_T const & dof, REAL* forceAtom_in, REAL* forceAtom_out, unsigned const idx_protein,
		typename std::enable_if< std::is_same<DOF_T,DOF_6D_Modes<REAL>>::value, void>::type* dummy = 0)
		//typename std::enable_if< std::is_same<DOF_T,DOF_6D_Modes<REAL>>, void>::type* dummy1 = 0)
		{
	if(PROTEIN_T == 0){
			forceAtom_out[0] = forceAtom_in[0];
			forceAtom_out[1] = forceAtom_in[1];
			forceAtom_out[2] = forceAtom_in[2];
		}
	else{
		RotMat<REAL> const rotMatInv = euler2rotmat(dof._6D.ang.x, dof._6D.ang.y, dof._6D.ang.z).getInv();
		forceAtom_out[0] = rotMatInv[0] * forceAtom_in[0] + rotMatInv[1] * forceAtom_in[1] + rotMatInv[2] * forceAtom_in[2];
		forceAtom_out[1] = rotMatInv[3] * forceAtom_in[0] + rotMatInv[4] * forceAtom_in[1] + rotMatInv[5] * forceAtom_in[2];
		forceAtom_out[2] = rotMatInv[6] * forceAtom_in[0] + rotMatInv[7] * forceAtom_in[1] + rotMatInv[8] * forceAtom_in[2];
	}
}

template <typename REAL, typename DOF_T, int PROTEIN_T>
__device__ inline void rotateForce( DOF_T const & dof,  REAL* forceAtom_in, REAL* forceAtom_out, unsigned const idx_protein,
		typename std::enable_if<  std::is_same<DOF_T,DOF_6D<REAL>>::value, void>::type* dummy = 0 ){
	if(PROTEIN_T == 0){
		forceAtom_out[0] = forceAtom_in[0];
		forceAtom_out[1] = forceAtom_in[1];
		forceAtom_out[2] = forceAtom_in[2];
	}
}


template <typename REAL, typename DOF_T, unsigned blockSize>
__device__ inline void reduceModes( DOF_T const & dof,  REAL const* forceAtom, REAL const * xModes, REAL const * yModes, REAL const * zModes, unsigned const numModes, unsigned const idx_atom, REAL sum_mode[5],
		typename std::enable_if<std::is_same< DOF_T, DOF_MB_Modes<REAL> >::value, void>::type* dummy = 0 ){
	unsigned const idx_mode = (idx_atom + blockSize)* numModes;
	for ( unsigned mode = 0; mode < numModes; mode++){
		sum_mode[mode] -=   forceAtom[0]*xModes[idx_mode + mode]
						  + forceAtom[1]*yModes[idx_mode + mode]
						  + forceAtom[2]*zModes[idx_mode + mode];
	}
}

template <typename REAL, typename DOF_T, unsigned blockSize>
__device__ inline void reduceModes( DOF_T const & dof,  REAL const* forceAtom, REAL const * xModes, REAL const * yModes, REAL const * zModes, unsigned const numModes, unsigned const idx_atom, REAL sum_mode[5],
		typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0 ){
	unsigned const idx_mode = (idx_atom + blockSize)* numModes;
	for ( unsigned mode = 0; mode < numModes; mode++){
		sum_mode[mode] -=   forceAtom[0]*xModes[idx_mode + mode]
						  + forceAtom[1]*yModes[idx_mode + mode]
						  + forceAtom[2]*zModes[idx_mode + mode];
	}
}

template <typename REAL, typename DOF_T, unsigned blockSize>
__device__ inline void reduceModes( DOF_T const & dof,  REAL const* forceAtom, REAL const * xModes, REAL const * yModes, REAL const * zModes, unsigned const numModes, unsigned const idx_atom, REAL sum_mode[5],
		typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0 ){}


template <typename T,typename DOF_T, unsigned int blockSize, bool nIsPow2, int PROTEINTYPE >
__global__ void
blockReduceMode(unsigned const numAtoms, unsigned const numModes,  DOF_T* dofs,unsigned const idx_protein,
			T* xPos, T* yPos, T* zPos,
			T *xModes,T *yModes,T *zModes,
			T *d_fx, T *d_fy, T *d_fz,
			T *d_E,
			T *g_odata)

{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = threadIdx.x; // half the number of blocks
    unsigned int base = blockIdx.x * numAtoms;
    constexpr bool MODES = std::is_same<DOF_T,DOF_6D_Modes<T>>::value || std::is_same<DOF_T, DOF_MB_Modes<T>>::value;

    T sum_fx = 0;
    T sum_fy = 0;
    T sum_fz = 0;
    T sum_E = 0;
    T sum_torque[9] = {0};
    T sum_mode[5] = {0};


    while (i < numAtoms)
    {
    	unsigned int DOFidx = (base + i) / numAtoms;
    	auto const  dof = dofs[DOFidx];


		T f[3], x, y, z,f_mode[3];

		f[0] = d_fx[base + i];
		f[1] = d_fy[base + i];
		f[2] = d_fz[base + i];
		//
		x = xPos[base + i];
		y = yPos[base + i];
		z = zPos[base + i];
		sum_fx += f[0];
		sum_fy += f[1];
		sum_fz += f[2];
		sum_E += d_E[base + i];
		sum_torque[0] += x * f[0];
		sum_torque[1] += y * f[0];
		sum_torque[2] += z * f[0];
		sum_torque[3] += x * f[1];
		sum_torque[4] += y * f[1];
		sum_torque[5] += z * f[1];
		sum_torque[6] += x * f[2];
		sum_torque[7] += y * f[2];
		sum_torque[8] += z * f[2];

		rotateForce<T,DOF_T,PROTEINTYPE>(  dof, f, f_mode, idx_protein);
		reduceModes<T,DOF_T, 0>(  dof,  f_mode,  xModes, yModes,  zModes,  numModes, i,  sum_mode);

        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < numAtoms) {
			f[0] = d_fx[base + i + blockSize];
			f[1] = d_fy[base + i + blockSize];
			f[2] = d_fz[base + i + blockSize];
			x = xPos[i + blockSize];
			y = yPos[i + blockSize];
			z = zPos[i + blockSize];
			sum_fx += f[0];
			sum_fy += f[1];
			sum_fz += f[2];
			sum_E += d_E[base + i + blockSize];
			sum_torque[0] += x * f[0];
			sum_torque[1] += y * f[0];
			sum_torque[2] += z * f[0];
			sum_torque[3] += x * f[1];
			sum_torque[4] += y * f[1];
			sum_torque[5] += z * f[1];
			sum_torque[6] += x * f[2];
			sum_torque[7] += y * f[2];
			sum_torque[8] += z * f[2];


			rotateForce<T,DOF_T,PROTEINTYPE>(  dof, f, f_mode, idx_protein);
			reduceModes<T,DOF_T, blockSize>(  dof,  f_mode,  xModes, yModes,  zModes,  numModes, i,  sum_mode);

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
	if(MODES){
		for(int mode=0; mode < numModes; mode++){
			sdata[tid + (13+mode) * blockSize] =  sum_mode[mode];
		}
	}



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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13+mode)* blockSize] 	   = sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13+mode)* blockSize + 512];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13 + mode)* blockSize + 256];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13 + mode)* blockSize + 128];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13 + mode)* blockSize + 64];
			}
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
			if(MODES){
				for(int mode=0;mode<numModes;mode++){
					sum_mode[mode] += sdata[tid
										   + (13 + mode)* blockSize + 32];
				}
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
			if(MODES){
				for(int mode=0;mode<numModes;mode++){
					sum_mode[mode] +=__shfl_down(sum_mode[mode], offset);
				}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13 + mode)* blockSize + 32];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13+mode)* blockSize + 16];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13+mode)* blockSize + 8];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13 + mode)* blockSize + 4];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13 + mode)* blockSize + 2];
			}
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
		if(MODES){
			for(int mode=0;mode<numModes;mode++){
				sdata[tid + (13 + mode)* blockSize] =  sum_mode[mode] = sum_mode[mode]
																	   + sdata[tid + (13+mode)* blockSize + 1];
			}
		}

    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) {
		if(MODES){
			g_odata[0  + blockIdx.x*(13 + numModes)] = sum_fx;
			g_odata[1  + blockIdx.x*(13 + numModes)] = sum_fy;
			g_odata[2  + blockIdx.x*(13 + numModes)] = sum_fz;
			g_odata[3  + blockIdx.x*(13 + numModes)] = sum_E;
			g_odata[4  + blockIdx.x*(13 + numModes)] = sum_torque[0];
			g_odata[5  + blockIdx.x*(13 + numModes)] = sum_torque[1];
			g_odata[6  + blockIdx.x*(13 + numModes)] = sum_torque[2];
			g_odata[7  + blockIdx.x*(13 + numModes)] = sum_torque[3];
			g_odata[8  + blockIdx.x*(13 + numModes)] = sum_torque[4];
			g_odata[9  + blockIdx.x*(13 + numModes)] = sum_torque[5];
			g_odata[10 + blockIdx.x*(13 + numModes)] = sum_torque[6];
			g_odata[11 + blockIdx.x*(13 + numModes)] = sum_torque[7];
			g_odata[12 + blockIdx.x*(13 + numModes)] = sum_torque[8];

			for(int mode=0;mode<numModes;mode++){
				g_odata[13 + mode + blockIdx.x*(13 + numModes )] =  sum_mode[mode];
			}
		}
		else{
			g_odata[0  + blockIdx.x* 13] = sum_fx;
			g_odata[1  + blockIdx.x* 13] = sum_fy;
			g_odata[2  + blockIdx.x* 13] = sum_fz;
			g_odata[3  + blockIdx.x* 13] = sum_E;
			g_odata[4  + blockIdx.x* 13] = sum_torque[0];
			g_odata[5  + blockIdx.x* 13] = sum_torque[1];
			g_odata[6  + blockIdx.x* 13] = sum_torque[2];
			g_odata[7  + blockIdx.x* 13] = sum_torque[3];
			g_odata[8  + blockIdx.x* 13] = sum_torque[4];
			g_odata[9  + blockIdx.x* 13] = sum_torque[5];
			g_odata[10 + blockIdx.x* 13] = sum_torque[6];
			g_odata[11 + blockIdx.x* 13] = sum_torque[7];
			g_odata[12 + blockIdx.x* 13] = sum_torque[8];
		}

    }
}



template <typename T, typename DOF_T, int PROTEINTYPE>
void d_reduce(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtoms,
		const unsigned& numModes,
		const unsigned idx_protein,
		 DOF_T* dofs,
		T* xPos, T* yPos, T* zPos,
		T *xModes,T *yModes,T *zModes,
		T *d_fx, T *d_fy, T *d_fz,
		T *d_E,
		T *g_odata,
		const cudaStream_t& stream)

{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);
	int smemSize;
	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	if(std::is_same<DOF_T,DOF_6D_Modes<T>>::value ||std::is_same<DOF_T,DOF_MB_Modes<T>>::value ){
		smemSize = (threads <= 32) ? 2 * (13 + numModes) * threads * sizeof(T) : (13 + numModes) * threads * sizeof(T);
	}

	else{
		smemSize = (threads <= 32) ? 2 * 13 * threads * sizeof(T) : 13 * threads * sizeof(T);
	}
	// choose which of the optimized versions of reduction to launch



//	if (isPow2(max(numAtomsRec, numAtomsLig)))
		if (isPow2( numAtoms))
		{
			switch (threads)
			{
				case 1024:
					blockReduceMode<T, DOF_T, 1024,true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein, xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 512:
					blockReduceMode<T, DOF_T, 512, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs,idx_protein, xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 256:
					blockReduceMode<T, DOF_T, 256, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs,idx_protein, xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 128:
					blockReduceMode<T, DOF_T, 128, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein, xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 64:
					blockReduceMode<T, DOF_T,  64, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs,idx_protein, xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 32:
					blockReduceMode<T, DOF_T,  32, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 16:
					blockReduceMode<T, DOF_T,  16, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  8:
					blockReduceMode<T, DOF_T,   8, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  4:
					blockReduceMode<T, DOF_T,   4, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  2:
					blockReduceMode<T, DOF_T,   2, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  1:
					blockReduceMode<T, DOF_T,   1, true, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;
			}
		}
		else
		{
			switch (threads)
			{
				case 1024:
					blockReduceMode<T, DOF_T, 1024,false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 512:
					blockReduceMode<T, DOF_T, 512, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 256:
					blockReduceMode<T, DOF_T, 256, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 128:
					blockReduceMode<T, DOF_T, 128, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 64:
					blockReduceMode<T, DOF_T,  64, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 32:
					blockReduceMode<T, DOF_T,  32, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case 16:
					blockReduceMode<T, DOF_T,  16, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  8:
					blockReduceMode<T, DOF_T,   8, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes, d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  4:
					blockReduceMode<T, DOF_T,   4, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes,d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  2:
					blockReduceMode<T, DOF_T,   2, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes,d_fx, d_fy, d_fz, d_E, g_odata);
					break;

				case  1:
					blockReduceMode<T, DOF_T,   1, false, PROTEINTYPE><<< dimGrid, dimBlock, smemSize, stream >>>(numAtoms, numModes, dofs, idx_protein,xPos, yPos, zPos, xModes, yModes, zModes,d_fx, d_fy, d_fz, d_E, g_odata);
					break;
			}
		}

}


template<typename REAL, typename DOF_T, int PROTEINTYPE>
void deviceReduce(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<REAL>* protein,
		DOF_T* dofs,
		REAL* xPos, REAL* yPos, REAL* zPos,
		REAL *d_fx, REAL *d_fy, REAL *d_fz,
		REAL *d_E,
		REAL *d_out,
		const cudaStream_t& stream)
{
	/* we need at least twice the block size number of threads */
	const unsigned size = protein->numAtoms;
	const unsigned threads = (size < blockSize*2) ? nextPow2((size + 1)/ 2) : blockSize;
	/* each structure is reduced by one thread block */
	const unsigned blocks = numDOFs;
	d_reduce<REAL, DOF_T, PROTEINTYPE>(
		threads,
		blocks,
		size,
		protein->numModes,
		idx_protein,
		dofs,
		xPos,  yPos,  zPos,
		protein->xModes, protein->yModes, protein->zModes,
		d_fx,  d_fy,  d_fz,
		d_E,
		d_out,
		stream);
}

template
void deviceReduce< float, DOF_6D_Modes<float>, 0>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		float* xPos, float* yPos, float* zPos,
		float *d_fx, float *d_fy, float *d_fz,
		float *d_E,
		float *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< double, DOF_6D_Modes<double>, 0>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		double* xPos, double* yPos, double* zPos,
		double *d_fx, double *d_fy, double *d_fz,
		double *d_E,
		double *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< float, DOF_6D_Modes<float>, 1>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<float>* protein,
		DOF_6D_Modes<float>* dofs,
		float* xPos, float* yPos, float* zPos,
		float *d_fx, float *d_fy, float *d_fz,
		float *d_E,
		float *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< double, DOF_6D_Modes<double>, 1>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<double>* protein,
		DOF_6D_Modes<double>* dofs,
		double* xPos, double* yPos, double* zPos,
		double *d_fx, double *d_fy, double *d_fz,
		double *d_E,
		double *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< float, DOF_MB_Modes<float>, 1>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<float>* protein,
		DOF_MB_Modes<float>* dofs,
		float* xPos, float* yPos, float* zPos,
		float *d_fx, float *d_fy, float *d_fz,
		float *d_E,
		float *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< double, DOF_MB_Modes<double>, 1>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<double>* protein,
		DOF_MB_Modes<double>* dofs,
		double* xPos, double* yPos, double* zPos,
		double *d_fx, double *d_fy, double *d_fz,
		double *d_E,
		double *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< float, DOF_MB_Modes<float>, 0>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<float>* protein,
		DOF_MB_Modes<float>* dofs,
		float* xPos, float* yPos, float* zPos,
		float *d_fx, float *d_fy, float *d_fz,
		float *d_E,
		float *d_out,
		const cudaStream_t& stream);

template
void deviceReduce< double, DOF_MB_Modes<double>, 0>(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned idx_protein,
		d_Protein<double>* protein,
		DOF_MB_Modes<double>* dofs,
		double* xPos, double* yPos, double* zPos,
		double *d_fx, double *d_fy, double *d_fz,
		double *d_E,
		double *d_out,
		const cudaStream_t& stream);


}  // namespace as

