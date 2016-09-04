/*
 * interpolation.cu
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include "nativeTypesWrapper.h"
#include "DeviceIntrplGrid.h"
#include "DeviceProtein.h"
#include "macros.h"

namespace as {

template<typename REAL>
using d_IntrlpGrid = typename DeviceIntrplGrid<REAL>::Desc;

template<typename REAL>
using d_Protein = typename DeviceProtein<REAL>::Desc;

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

		data_out_x[idx] = static_cast<REAL>(pot.x);
		data_out_y[idx] = static_cast<REAL>(pot.y);
		data_out_z[idx] = static_cast<REAL>(pot.z);
		data_out_E[idx] = static_cast<REAL>(pot.w);
	}
}

template<typename REAL>
__global__ void d_outerPotForce(
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
		if (type != 0) {

			REAL x = data_in_x[idx];
			REAL y = data_in_y[idx];
			REAL z = data_in_z[idx];

			if (      ((x < grid.minDim.x || x > grid.maxDim.x)
					|| (y < grid.minDim.y || y > grid.maxDim.y)
					|| (z < grid.minDim.z || z > grid.maxDim.z))
					&&
					  ((x >= grid.minDim.x && x <= grid.maxDim.x)
					&& (y >= grid.minDim.y && y <= grid.maxDim.y)
					&& (z >= grid.minDim.z && z <= grid.maxDim.z)))
			{
				x = (x - grid.minDim.x) * grid.dVox_inv + 0.5f;
				y = (y - grid.minDim.y) * grid.dVox_inv + 0.5f;
				z = (z - grid.minDim.z) * grid.dVox_inv + 0.5f;

				float4 pot = tex3D<float4>(grid.texArrayLin[type], x, y, z); /** Interpolated value */

				REAL charge = prot.charge[idx % numAtoms];
				if (fabs(charge) > 0.001f) {
					float4 V_el = tex3D<float4>(grid.texArrayLin[0], x, y, z); /** Interpolated value */
					pot = pot + V_el * charge;
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


