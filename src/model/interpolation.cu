/*
 * interpolation.cu
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include "interpolation.h"



namespace as {

template <typename REAL>
__forceinline__ __device__ void getVoxelDevice(const d_IntrlpGrid<REAL>& grid, const unsigned &type,
		const float &x, const float &y,	const float &z,  VoxelOctet<float>& voxelOct, unsigned const i)
{
	unsigned idxX = (unsigned) floor(
			(x - grid.minDim.x) * grid.dVox_inv);
	unsigned idxY = (unsigned) floor(
			(y - grid.minDim.y) * grid.dVox_inv);
	unsigned idxZ = (unsigned) floor(
			(z - grid.minDim.z) * grid.dVox_inv);

	// compute absolute position of the vertices
	voxelOct.min.x = idxX * grid.dVox + grid.minDim.x;
	voxelOct.min.y = idxY * grid.dVox + grid.minDim.y;
	voxelOct.min.z = idxZ * grid.dVox + grid.minDim.z;
	voxelOct.max.x = voxelOct.min.x + grid.dVox;
	voxelOct.max.y = voxelOct.min.y + grid.dVox;
	voxelOct.max.z = voxelOct.min.z + grid.dVox;

	/** use of non-normalized coordinates */

	float idxNx =  0.5 +(float) idxX;  //
	float idxNy =  0.5 +(float) idxY;  //
	float idxNz =  0.5 +(float) idxZ;  //
	//printf("%d %f %f %f\n",i,idxNx,idxNy,idxNz);
	//printf("%d \t %0.1f %0.1f %0.1f \t %f %f %f %f\n",i,idxNx,idxNy,idxNz,voxelOct.data[0][0][0].x,voxelOct.data[0][0][0].y,voxelOct.data[0][0][0].z,voxelOct.data[0][0][0].w );
	voxelOct.data[0][0][0] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy, idxNz);
	voxelOct.data[1][0][0] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy, idxNz);
	voxelOct.data[0][1][0] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy + 1, idxNz);
	voxelOct.data[1][1][0] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy + 1, idxNz);

	voxelOct.data[0][0][1] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy, idxNz + 1);
	voxelOct.data[1][0][1] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy, idxNz + 1);
	voxelOct.data[0][1][1] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy + 1, idxNz + 1);
	voxelOct.data[1][1][1] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy + 1, idxNz + 1);

}


/*
** @brief: function body for a trilinear interpolation.
*/

__forceinline__ __host__ __device__ void trilinearInterpolation(const float &x,
		const float &y, const float &z, const VoxelOctet<float> &voxelOct,
		const float &voxelVol_inv, float4 &V)
{
	/* for operator overloading of *,+,-,/ for cuda types (REAL4)
	 * they are defined in asUtils/cudaMath*/

	 //printf("%.20f %.20f %.20f\n",x,y,z);
	float3 pos = make_float3(x, y, z);
	float3 pos_m_posMin ;
	pos_m_posMin.x = pos.x - voxelOct.min.x;
	pos_m_posMin.y = pos.y - voxelOct.min.y;
	pos_m_posMin.z = pos.z - voxelOct.min.z;
	float3 posMax_m_pos ;
	posMax_m_pos.x = voxelOct.max.x - pos.x;
	posMax_m_pos.y = voxelOct.max.y - pos.y;
	posMax_m_pos.z = voxelOct.max.z - pos.z;

	float tmpMax_xy = (posMax_m_pos.x) * (posMax_m_pos.y);
	float tmpMin_xy = (pos_m_posMin.x) * (pos_m_posMin.y);

	V = voxelOct.data[0][0][0] * (tmpMax_xy * (posMax_m_pos.z))
			+ voxelOct.data[1][0][0]
					* ((pos_m_posMin.x) * (posMax_m_pos.y) * (posMax_m_pos.z))
			+ voxelOct.data[0][1][0]
					* ((posMax_m_pos.x) * (pos_m_posMin.y) * (posMax_m_pos.z))
			+ voxelOct.data[0][0][1] * (tmpMax_xy * (pos_m_posMin.z))
			+ voxelOct.data[1][0][1]
					* ((pos_m_posMin.x) * (posMax_m_pos.y) * (pos_m_posMin.z))
			+ voxelOct.data[0][1][1]
					* ((posMax_m_pos.x) * (pos_m_posMin.y) * (pos_m_posMin.z))
			+ voxelOct.data[1][1][0] * (tmpMin_xy * (posMax_m_pos.z))
			+ voxelOct.data[1][1][1] * (tmpMin_xy * (pos_m_posMin.z));

	V = V * voxelVol_inv;
	return;
}
template<typename REAL>
__forceinline__ __device__ float4 Intrpl3D(const d_IntrlpGrid<REAL>& grid, const unsigned& type, const float &x, const float &y,
		const float &z, unsigned const i)
{

	VoxelOctet<float> voxelOct;
	getVoxelDevice<REAL>(grid, type, x, y, z, voxelOct,i);
	float4 V;
	trilinearInterpolation(x, y, z, voxelOct, grid.voxelVol_inv, V);
	return V;
}

template<typename REAL>
__global__ void d_innerPotForce (
		const d_IntrlpGrid<REAL> grid,
		const d_Protein<REAL> prot,
		const unsigned numDOFs,
		const REAL* data_in_x,
		const REAL* data_in_y,
		const REAL* data_in_z,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{
	using real4_t = typename TypeWrapper<REAL>::real4_t;
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned numAtoms = prot.numAtoms;
	if (idx < numAtoms*numDOFs) {
		unsigned type = prot.mappedType[idx % numAtoms];
		float4 pot = {0,0,0,0};
		if (type != 0) {

			REAL x = data_in_x[idx];
			REAL y = data_in_y[idx];
			REAL z = data_in_z[idx];
			bool test =false;

			if ((x >= grid.minDim.x && x <= grid.maxDim.x)
					&& (y >= grid.minDim.y && y <= grid.maxDim.y)
					&& (z >= grid.minDim.z && z <= grid.maxDim.z))
			{
				if(test){

					x = (x - grid.minDim.x) * grid.dVox_inv + 0.5;
					y = (y - grid.minDim.y) * grid.dVox_inv + 0.5;
					z = (z - grid.minDim.z) * grid.dVox_inv + 0.5;

					pot = tex3D<float4>(grid.texArrayLin[type], x, y, z); /** Interpolated value */

					REAL charge = prot.charge[idx % numAtoms];
					if (fabs(charge) > 0.001f) {
						float4 V_el = tex3D<float4>(grid.texArrayLin[0], x, y, z); /** Interpolated value */
						pot = pot + V_el * charge;
					}
				}else{

					//pot = interpolate2<REAL>(  grid,type,  x, y,  z,idx);
					pot = Intrpl3D<REAL>(grid, type, (float)x,  (float)y,  (float)z,idx);

					//printf("%f %f\n", pot.y,test.y);
					//pot = tex3Dfetch<float4>(grid.texArrayLin[type], {idxX, idxY, idxZ,1}); /** Interpolated value */

					REAL charge = prot.charge[idx % numAtoms];
					if (fabs(charge) > 0.001f) {

						float4 V_el = Intrpl3D<REAL>(grid, 0, x, y, z,idx);
						//float4 V_el= interpolate2<REAL>(  grid,0,  x, y,  z,idx);
						pot = pot + V_el * charge;
					}
				}
			}
		}

		data_out_x[idx] = static_cast<REAL>(pot.x);
		data_out_y[idx] = static_cast<REAL>(pot.y);
		data_out_z[idx] = static_cast<REAL>(pot.z);
		data_out_E[idx] = static_cast<REAL>(pot.w);
	}
}



template<typename REAL>
__global__ void d_outerPotForce(
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		const d_Protein<REAL> prot,
		const unsigned numDOFs,
		const REAL* data_in_x,
		const REAL* data_in_y,
		const REAL* data_in_z,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{

	using real4_t = typename TypeWrapper<REAL>::real4_t;
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

	//DEBUG
//	if (idx == 0) {
//		printf("%f %f %f %f %f %f\n" ,
//				grid.minDim.x, grid.minDim.y, grid.minDim.z,
//				grid.maxDim.x, grid.maxDim.y, grid.maxDim.z);
//	}


	const unsigned numAtoms = prot.numAtoms;
	if (idx < numAtoms*numDOFs) {
		unsigned type = prot.mappedType[idx % numAtoms];
		if (type != 0) {

			REAL x = data_in_x[idx];
			REAL y = data_in_y[idx];
			REAL z = data_in_z[idx];

//			if (idx < 50) {
//				printf("%u %f %f %f %u\n" ,
//						idx, x, y, z, type);
//			}

			if (      ((x < inner.minDim.x || x > inner.maxDim.x)
					|| (y < inner.minDim.y || y > inner.maxDim.y)
					|| (z < inner.minDim.z || z > inner.maxDim.z))
					&&
					  ((x >= outer.minDim.x && x <= outer.maxDim.x)
					&& (y >= outer.minDim.y && y <= outer.maxDim.y)
					&& (z >= outer.minDim.z && z <= outer.maxDim.z)))
			{
				bool test =false;
				float4 pot{0,0,0,0};
				if(test){
				x = (x - outer.minDim.x) * outer.dVox_inv + 0.5f;
				y = (y - outer.minDim.y) * outer.dVox_inv + 0.5f;
				z = (z - outer.minDim.z) * outer.dVox_inv + 0.5f;

				pot = tex3D<float4>(outer.texArrayLin[type], x, y, z); /** Interpolated value */

				REAL charge = prot.charge[idx % numAtoms];
				if (fabs(charge) > 0.001f) {
					float4 V_el = tex3D<float4>(outer.texArrayLin[0], x, y, z); /** Interpolated value */
					pot = pot + V_el * charge;
				}

//				if (idx < 20) {
//					printf("%u %f %f %f %f %f %f\n" ,
//							idx, pot.x, pot.y, pot.z, pot.w);
//				}
				}else{


					pot = Intrpl3D<REAL>(outer, type, (float)x,  (float)y,  (float)z,idx);
					//pot = interpolate2<REAL>(  outer, type,  x, y,  z,idx);


					REAL charge = prot.charge[idx % numAtoms];
					if (fabs(charge) > 0.001f) {

						float4 V_el = Intrpl3D<REAL>(outer, 0, x, y, z,idx);
						//float4 V_el = interpolate2<REAL>(  outer, 0,  x, y,  z,idx);
						pot = pot + V_el * charge;
					}
				}

				data_out_x[idx] = static_cast<REAL>(pot.x);
				data_out_y[idx] = static_cast<REAL>(pot.y);
				data_out_z[idx] = static_cast<REAL>(pot.z);
				data_out_E[idx] = static_cast<REAL>(pot.w);
			}
		}
	}
}

template<typename REAL>
void d_potForce (
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<REAL>& inner,
		const d_IntrlpGrid<REAL>& outer,
		const d_Protein<REAL>& prot,
		const unsigned& numDOFs,
		const REAL* data_in_x,
		const REAL* data_in_y,
		const REAL* data_in_z,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{
	cudaVerifyKernel((
			d_innerPotForce<<<gridSize, blockSize, 0, stream>>> (
				inner,
				prot,
				numDOFs,
				data_in_x,
				data_in_y,
				data_in_z,
				data_out_x,
				data_out_y,
				data_out_z,
				data_out_E
			)
		));

	cudaVerifyKernel((
			d_outerPotForce<<<gridSize, blockSize, 0, stream>>> (
				inner,
				outer,
				prot,
				numDOFs,
				data_in_x,
				data_in_y,
				data_in_z,
				data_out_x,
				data_out_y,
				data_out_z,
				data_out_E
			)
		));
}

template
void d_potForce<float> (
	unsigned blockSize,	unsigned gridSize, const cudaStream_t &stream,
	const d_IntrlpGrid<float>& inner, const d_IntrlpGrid<float>& outer, const d_Protein<float>& prot,
	const unsigned& numDOFs,
	const float* data_in_x, const float* data_in_y, const float* data_in_z,
	float* data_out_x, float* data_out_y, float* data_out_z, float* data_out_E);

template
void d_potForce<double> (
	unsigned blockSize,	unsigned gridSize, const cudaStream_t &stream,
	const d_IntrlpGrid<double>& inner, const d_IntrlpGrid<double>& outer, const d_Protein<double>& prot,
	const unsigned& numDOFs,
	const double* data_in_x, const double* data_in_y, const double* data_in_z,
	double* data_out_x, double* data_out_y, double* data_out_z, double* data_out_E);

}  // namespace as


