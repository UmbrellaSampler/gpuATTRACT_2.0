/*
 * interpolation_modes.cu
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include "interpolation_MB_modes.h"

#include "nativeTypesWrapper.h"
#include "DeviceIntrplGrid.h"
#include "DeviceProtein.h"
#include "macros.h"

namespace as {

template<typename REAL>
__global__ void d_innerPotForce (
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned int ligIdx,
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
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		unsigned type = prot.mappedType[idx % numAtoms];
		float4 pot = {0,0,0,0};
		if (type != 0) {

			REAL x = data_in_x[idx];
			REAL y = data_in_y[idx];
			REAL z = data_in_z[idx];

			if ((x >= grid.minDim.x && x <= grid.maxDim.x)
					&& (y >= grid.minDim.y && y <= grid.maxDim.y)
					&& (z >= grid.minDim.z && z <= grid.maxDim.z))
			{
				x = (x - grid.minDim.x) * grid.dVox_inv + 0.5;
				y = (y - grid.minDim.y) * grid.dVox_inv + 0.5;
				z = (z - grid.minDim.z) * grid.dVox_inv + 0.5;

				pot = tex3D<float4>(grid.texArrayLin[type], x, y, z); /** Interpolated value */

				REAL charge = prot.charge[idx % numAtoms];
				if (fabs(charge) > 0.001f) {
					float4 V_el = tex3D<float4>(grid.texArrayLin[0], x, y, z); /** Interpolated value */
					pot = pot + V_el * charge;
				}
			}
		}

		Vec3<REAL> force(pot.x, pot.y, pot.z);
		const RotMat<REAL> rotMat = euler2rotmat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);
		force = rotMat * force;

		data_out_x[idx] += static_cast<REAL>(force.x);
		data_out_y[idx] += static_cast<REAL>(force.y);
		data_out_z[idx] += static_cast<REAL>(force.z);
		data_out_E[idx] += static_cast<REAL>(pot.w);
	}
}

template<typename REAL>
__global__ void d_outerPotForce(
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned int ligIdx,
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
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
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

				x = (x - outer.minDim.x) * outer.dVox_inv + 0.5f;
				y = (y - outer.minDim.y) * outer.dVox_inv + 0.5f;
				z = (z - outer.minDim.z) * outer.dVox_inv + 0.5f;

				float4 pot = tex3D<float4>(outer.texArrayLin[type], x, y, z); /** Interpolated value */

				REAL charge = prot.charge[idx % numAtoms];
				if (fabs(charge) > 0.001f) {
					float4 V_el = tex3D<float4>(outer.texArrayLin[0], x, y, z); /** Interpolated value */
					pot = pot + V_el * charge;
				}

//				if (idx < 20) {
//					printf("%u %f %f %f %f %f %f\n" ,
//							idx, pot.x, pot.y, pot.z, pot.w);
//				}

				Vec3<REAL> force(pot.x, pot.y, pot.z);
				const RotMat<REAL> rotMat = euler2rotmat(dof._6D[ligIdx].ang.x, dof._6D[ligIdx].ang.y, dof._6D[ligIdx].ang.z);
				force = rotMat * force;

				data_out_x[idx] += static_cast<REAL>(force.x);
				data_out_y[idx] += static_cast<REAL>(force.y);
				data_out_z[idx] += static_cast<REAL>(force.z);
				data_out_E[idx] += static_cast<REAL>(pot.w);
			}
		}
	}
}

template<typename REAL>
void d_potForce (
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		DOF_6D_MB_Modes<REAL>* dofs,
		unsigned int ligIdx,
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
				dofs,
				ligIdx,
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
				dofs,
				ligIdx,
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
	DOF_6D_MB_Modes<float>* dofs,
	unsigned int ligIdx,
	const d_IntrlpGrid<float>& inner, const d_IntrlpGrid<float>& outer, const d_Protein<float>& prot,
	const unsigned& numDOFs,
	const float* data_in_x, const float* data_in_y, const float* data_in_z,
	float* data_out_x, float* data_out_y, float* data_out_z, float* data_out_E);

template
void d_potForce<double> (
	unsigned blockSize,	unsigned gridSize, const cudaStream_t &stream,
	DOF_6D_MB_Modes<double>* dofs,
	unsigned int ligIdx,
	const d_IntrlpGrid<double>& inner, const d_IntrlpGrid<double>& outer, const d_Protein<double>& prot,
	const unsigned& numDOFs,
	const double* data_in_x, const double* data_in_y, const double* data_in_z,
	double* data_out_x, double* data_out_y, double* data_out_z, double* data_out_E);

}  // namespace as


