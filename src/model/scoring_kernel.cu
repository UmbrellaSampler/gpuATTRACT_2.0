#include "scoring_kernel.h"
#include "Protein.h"
#include "IntrplGrid.h"
#include "nativeTypesFunctions.h"
#include "nativeTypesMath.h"
#include "VoxelOctet.h"
#include "trilinIntrpl.h"
#include <cmath>

#ifdef CUDA
#include "nativeTypesWrapper.h"
#include "DeviceIntrplGrid.h"
#include "DeviceProtein.h"
#include "macros.h"

#endif
namespace as{
template <typename T>
__host__ __device__
__forceinline__ T lerp(T v0, T v1, T t) {
    return fma(t, v1, fma(-t, v0, v0));
}


__host__ __device__
__forceinline__ float4 lerp4f(float4 v0, float4 v1, float t) {
    return make_float4( lerp<float>(v0.x, v1.x, t), lerp<float>(v0.y, v1.y, t), lerp<float>(v0.z, v1.z, t), lerp<float>(v0.w, v1.w, t) );
}


template<typename REAL>
__host__ __device__
__forceinline__ float4 interpolate2( const d_IntrlpGrid<REAL>& grid,unsigned const type, REAL x, REAL y, REAL z, unsigned const i){
//TEMP
//	unsigned idxX1 = (unsigned) floor(
//				(x - grid.minDim.x) * grid.dVox_inv);
//		unsigned idxY1 = (unsigned) floor(
//				(y - grid.minDim.y) * grid.dVox_inv);
//		unsigned idxZ1 = (unsigned) floor(
//				(z - grid.minDim.z) * grid.dVox_inv);
//
//		// compute absolute position of the vertices
//		VoxelOctet<float> voxelOct;
//		voxelOct.min.x = idxX1 * grid.dVox + grid.minDim.x;
//		voxelOct.min.y = idxY1 * grid.dVox + grid.minDim.y;
//		voxelOct.min.z = idxZ1 * grid.dVox + grid.minDim.z;
//		voxelOct.max.x = voxelOct.min.x + grid.dVox;
//		voxelOct.max.y = voxelOct.min.y + grid.dVox;
//		voxelOct.max.z = voxelOct.min.z + grid.dVox;
//		voxelOct.data[0][0][0] = tex3D<float4>(grid.texArrayPt[type], idxX1, idxY1, idxZ1);
//		/** use of non-normalized coordinates */
//
//		float idxNx1 =  0.5 +(float) idxX1;  //
//		float idxNy1 =  0.5 +(float) idxY1;  //
//		float idxNz1 =  0.5 +(float) idxZ1;  //
//
//		float3 pos = make_float3(x, y, z);
//			float3 pos_m_posMin ;
//			pos_m_posMin.x = pos.x - voxelOct.min.x;
//			pos_m_posMin.y = pos.y - voxelOct.min.y;
//			pos_m_posMin.z = pos.z - voxelOct.min.z;
//			float3 posMax_m_pos ;
//			posMax_m_pos.x = voxelOct.max.x - pos.x;
//			posMax_m_pos.y = voxelOct.max.y - pos.y;
//			posMax_m_pos.z = voxelOct.max.z - pos.z;
//
//			float tmpMax_xy = (posMax_m_pos.x) * (posMax_m_pos.y);
//			float tmpMin_xy = (pos_m_posMin.x) * (pos_m_posMin.y);
//
//			float4 V = voxelOct.data[0][0][0] * (tmpMax_xy * (posMax_m_pos.z))
//					+ voxelOct.data[1][0][0]
//							* ((pos_m_posMin.x) * (posMax_m_pos.y) * (posMax_m_pos.z))
//					+ voxelOct.data[0][1][0]
//							* ((posMax_m_pos.x) * (pos_m_posMin.y) * (posMax_m_pos.z))
//					+ voxelOct.data[0][0][1] * (tmpMax_xy * (pos_m_posMin.z))
//					+ voxelOct.data[1][0][1]
//							* ((pos_m_posMin.x) * (posMax_m_pos.y) * (pos_m_posMin.z))
//					+ voxelOct.data[0][1][1]
//							* ((posMax_m_pos.x) * (pos_m_posMin.y) * (pos_m_posMin.z))
//					+ voxelOct.data[1][1][0] * (tmpMin_xy * (posMax_m_pos.z))
//					+ voxelOct.data[1][1][1] * (tmpMin_xy * (pos_m_posMin.z));



		//TEMP
	x = (x - grid.minDim.x) * grid.dVox_inv;
	y = (y - grid.minDim.y) * grid.dVox_inv;
	z = (z - grid.minDim.z) * grid.dVox_inv;

	unsigned const idxX = (unsigned) floor(x);
	unsigned const idxY = (unsigned) floor(y);
	unsigned const idxZ = (unsigned) floor(z);

	REAL const a = x - (REAL)idxX;
	REAL const b = y - (REAL)idxY;
	REAL const c = z - (REAL)idxZ;
	float4 data[2][2][2];
	data[0][0][0] = tex3D<float4>(grid.texArrayPt[type], idxX, idxY, idxZ);
	data[0][0][1] = tex3D<float4>(grid.texArrayPt[type], idxX, idxY, idxZ + 1);
	data[0][1][1] = tex3D<float4>(grid.texArrayPt[type], idxX, idxY + 1, idxZ + 1);
	data[0][1][0] = tex3D<float4>(grid.texArrayPt[type], idxX, idxY + 1, idxZ);
	data[1][1][0] = tex3D<float4>(grid.texArrayPt[type], idxX + 1, idxY + 1, idxZ);
	data[1][1][1] = tex3D<float4>(grid.texArrayPt[type], idxX + 1, idxY + 1, idxZ + 1);
	data[1][0][1] = tex3D<float4>(grid.texArrayPt[type], idxX + 1, idxY, idxZ + 1);
	data[1][0][0] = tex3D<float4>(grid.texArrayPt[type], idxX + 1, idxY, idxZ);

	float4 result =	lerp4f(
					lerp4f(
						lerp4f(data[0][0][0],data[0][0][1],c),
						lerp4f(data[0][1][0],data[0][1][1],c),
						b),
					lerp4f(
						lerp4f(data[1][0][0],data[1][0][1],c),
						lerp4f(data[1][1][0],data[1][1][1],c),
						b),
					a);


	//printf("%d %d %d %d     %f %f %f\n",i,idxX,idxY,idxZ,idxNx1,idxNy1,idxNz1);
	//printf("%d %f \n",i,data[0][0][0].x - voxelOct.data[0][0][0].x);
	//printf("%d %f %f\n",i,result.x, V.x);


return result;


}





template<typename REAL>
__device__ __forceinline__ void PotForce_device(
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		const d_Protein<REAL> prot,
		const unsigned numDOFs,
		const unsigned idx,
		const REAL x,
		const REAL y,
		const REAL z,
		float4 & data_out
		)
{

	using real4_t = typename TypeWrapper<REAL>::real4_t;


	const unsigned numAtoms = prot.numAtoms;
	unsigned type = prot.mappedType[idx % numAtoms];
	REAL charge = prot.charge[idx % numAtoms];
	if (type != 0) {
		if ((x >= inner.minDim.x && x <= inner.maxDim.x)
		 && (y >= inner.minDim.y && y <= inner.maxDim.y)
		 && (z >= inner.minDim.z && z <= inner.maxDim.z)){
			 gridForce( inner, x, y, z,idx, type,charge, data_out);
		}

		else if ( ((x >= outer.minDim.x && x <= outer.maxDim.x)
				&& (y >= outer.minDim.y && y <= outer.maxDim.y)
				&& (z >= outer.minDim.z && z <= outer.maxDim.z)))
		{
			gridForce( outer, x, y, z, idx,type,charge, data_out);
		}
	}
}


template<typename REAL>
__device__ __forceinline__ void gridForce(
		const d_IntrlpGrid<REAL> grid,
		 REAL x,
		 REAL y,
		 REAL z,
		const unsigned idx,
		const unsigned type,
		REAL charge,
		float4& data_out
		)
{

	using real4_t = typename TypeWrapper<REAL>::real4_t;

	bool test =false;
	if(test){
		x = (x - grid.minDim.x) * grid.dVox_inv + 0.5f;
		y = (y - grid.minDim.y) * grid.dVox_inv + 0.5f;
		z = (z - grid.minDim.z) * grid.dVox_inv + 0.5f;
		data_out = tex3D<float4>(grid.texArrayLin[type], x, y, z); /** Interpolated value */


		if (fabs(charge) > 0.001f) {
			float4 V_el = tex3D<float4>(grid.texArrayLin[0], x, y, z); /** Interpolated value */
			data_out = data_out + V_el * charge;
		}
	}else{

		data_out = interpolate2(  grid, type,  x, y,  z,idx);
		//REAL charge = prot.charge[idx % numAtoms];
		if (fabs(charge) > 0.001f) {


			float4 V_el = interpolate2(  grid, 0,  x, y,  z,idx);
			data_out = data_out + V_el * charge;
		}
	}
}

template<typename REAL, typename DOF_T >
__global__ void scoring_kernel(
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		d_Protein<REAL> const  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{
	using real4_t = typename TypeWrapper<REAL>::real4_t;
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned numAtoms = protein.numAtoms;
	if (idx < numAtoms*numDOFs) {
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		REAL x_trafo,y_trafo,z_trafo;
		float4 potForce;
		d_DOFPos_device( protein, dof, idx,	type_protein,
				 buffer_defoX[idx],  buffer_defoY[idx], buffer_defoZ[idx],
				 x_trafo,y_trafo,z_trafo
				);

		PotForce_device(inner, outer, protein, numDOFs, idx, x_trafo, y_trafo, z_trafo, potForce);
		 buffer_trafoX[idx] = x_trafo;
		 buffer_trafoY[idx] = y_trafo;
		 buffer_trafoZ[idx] = z_trafo;
	}
}


template<typename REAL, typename DOF_T>
 void d_score(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		d_Protein<REAL> const  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{
	cudaVerifyKernel((
			scoring_kernel<<<gridSize, blockSize, 0, stream>>> (
			inner,
			outer,
			protein,
			dofs,
			numDOFs,
			type_protein,
			buffer_defoX,
			buffer_defoY,
			buffer_defoZ,
			buffer_trafoX,
			buffer_trafoY,
			buffer_trafoZ,
			data_out_x,
			data_out_y,
			data_out_z,
			data_out_E))
		);
}
template
 void d_score<float, DOF_6D<float>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<float> inner,
		const d_IntrlpGrid<float> outer,
		d_Protein<float> const  protein,
		DOF_6D<float> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		float* buffer_defoX,
		float* buffer_defoY,
		float* buffer_defoZ,
		float* buffer_trafoX,
		float* buffer_trafoY,
		float* buffer_trafoZ,
		float* data_out_x,
		float* data_out_y,
		float* data_out_z,
		float* data_out_E);

template
 void d_score<double, DOF_6D<double>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<double> inner,
		const d_IntrlpGrid<double> outer,
		d_Protein<double> const  protein,
		DOF_6D<double> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		double* buffer_defoX,
		double* buffer_defoY,
		double* buffer_defoZ,
		double* buffer_trafoX,
		double* buffer_trafoY,
		double* buffer_trafoZ,
		double* data_out_x,
		double* data_out_y,
		double* data_out_z,
		double* data_out_E);

template
 void d_score<float, DOF_6D_Modes<float>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<float> inner,
		const d_IntrlpGrid<float> outer,
		d_Protein<float> const  protein,
		DOF_6D_Modes<float> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		float* buffer_defoX,
		float* buffer_defoY,
		float* buffer_defoZ,
		float* buffer_trafoX,
		float* buffer_trafoY,
		float* buffer_trafoZ,
		float* data_out_x,
		float* data_out_y,
		float* data_out_z,
		float* data_out_E);

template
 void d_score<double, DOF_6D_Modes<double>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<double> inner,
		const d_IntrlpGrid<double> outer,
		d_Protein<double> const  protein,
		DOF_6D_Modes<double> const * dofs,
		unsigned const numDOFs,
		unsigned const type_protein,
		double* buffer_defoX,
		double* buffer_defoY,
		double* buffer_defoZ,
		double* buffer_trafoX,
		double* buffer_trafoY,
		double* buffer_trafoZ,
		double* data_out_x,
		double* data_out_y,
		double* data_out_z,
		double* data_out_E);

}
