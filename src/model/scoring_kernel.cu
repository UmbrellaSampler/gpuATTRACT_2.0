#include "scoring_kernel.h"

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
			&& (z >= inner.minDim.z && z <= inner.maxDim.z))
		{
			 gridForce( inner, x, y, z,idx, type,charge, data_out);

		}

		if (      ((x < inner.minDim.x || x > inner.maxDim.x)
							|| (y < inner.minDim.y || y > inner.maxDim.y)
							|| (z < inner.minDim.z || z > inner.maxDim.z))
							&&
							  ((x >= outer.minDim.x && x <= outer.maxDim.x)
							&& (y >= outer.minDim.y && y <= outer.maxDim.y)
							&& (z >= outer.minDim.z && z <= outer.maxDim.z))){
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

	bool test =true;
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
		if (fabs(charge) > 0.001f) {
			float4 V_el = interpolate2(  grid, 0,  x, y,  z,idx);
			data_out = data_out + V_el * charge;
		}
	}
}


template<typename REAL, typename DOF_T >
__global__ void deform_kernel(
		d_Protein<REAL> const  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const idx_protein,
		REAL* posX,
		REAL* posY,
		REAL* posZ,

		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ
		)
{
	using real4_t = typename TypeWrapper<REAL>::real4_t;
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned numAtoms = protein.numAtoms;
	if (idx < numAtoms*numDOFs) {
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		const unsigned int num_atoms = protein.numAtoms;

		/* load DOF from global memory */
		unsigned atomIdx = idx % num_atoms;

		Vec3<REAL> posAtom(	posX[atomIdx],
							posY[atomIdx],
							posZ[atomIdx]);
		deform< REAL, DOF_T>(	 dof,posAtom, protein, atomIdx, idx_protein,
				buffer_defoX[idx], buffer_defoY[idx], buffer_defoZ[idx]);
	}
}



template<typename REAL, typename DOF_T >
__global__ void transform_kernel(
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const numAtoms,
		d_Protein<REAL> proteinCenter,
		d_Protein<REAL> proteinPartner,
		unsigned const idx_proteinCenter,
		unsigned const idx_proteinPartner,
//		Vec3<REAL> pivotCenter,
//		Vec3<REAL> pivotPartner,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ
		)
{
	using real4_t = typename TypeWrapper<REAL>::real4_t;
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < numAtoms*numDOFs) {
		unsigned DOFidx = idx / numAtoms;
		auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;
		Vec3<REAL> posAtomCenter(	buffer_defoX[idx],
									buffer_defoY[idx],
									buffer_defoZ[idx]);
		//Vec3<REAL> deformVec(0.0f);

//		deform< REAL, DOF_T>(	 dof,posAtomCenter, proteinCenter, atomIdx, idx_proteinCenter,
//				deformVec.x, deformVec.y, deformVec.z);


		device_transform(dof, posAtomCenter,proteinCenter.pivot, proteinPartner.pivot, idx_proteinCenter, idx_proteinPartner,
		       		 buffer_trafoX[idx],  buffer_trafoY[idx], buffer_trafoZ[idx]);

		//DEBUG
//		buffer_trafoX[idx] += pivotCenter.x;
//		buffer_trafoY[idx] += pivotCenter.y;
//		buffer_trafoZ[idx] += pivotCenter.z;

	}
}

template<typename REAL>
__global__ void transform_kernel_new(
		dofMB<REAL>* dofSupport,
		unsigned const numDOFs,
		unsigned const numAtoms,
		unsigned const numProtein,
		unsigned const idx_proteinCenter,
		unsigned const idx_proteinPartner,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ
		)
{
	using real4_t = typename TypeWrapper<REAL>::real4_t;
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx < numAtoms*numDOFs) {
		unsigned DOFidx = idx / numAtoms;
		//auto dof = dofs[DOFidx];
		unsigned atomIdx = idx % numAtoms;
		Vec3<REAL> posAtomCenter(	buffer_defoX[idx],
									buffer_defoY[idx],
									buffer_defoZ[idx]);
		//Vec3<REAL> deformVec(0.0f);

//		deform< REAL, DOF_T>(	 dof,posAtomCenter, proteinCenter, atomIdx, idx_proteinCenter,
//				deformVec.x, deformVec.y, deformVec.z);

		dofMB<REAL> dof = dofSupport[numDOFs * numProtein * idx_proteinCenter + numProtein *DOFidx + idx_proteinPartner ];

		posAtomCenter = dof.rotMat * posAtomCenter;
		posAtomCenter += dof.translation;
//		device_transform(dof, posAtomCenter,proteinCenter.pivot, proteinPartner.pivot, idx_proteinCenter, idx_proteinPartner,
//		       		 buffer_trafoX[idx],  buffer_trafoY[idx], buffer_trafoZ[idx]);

		//DEBUG
		buffer_trafoX[idx] = posAtomCenter.x;
		buffer_trafoY[idx] = posAtomCenter.y;
		buffer_trafoZ[idx] = posAtomCenter.z;

	}
}

//template<typename REAL, typename DOF_T >
//__global__ void transform_kernel(
//		DOF_T const * dofs,
//		unsigned const numDOFs,
//		unsigned const numAtoms,
//		d_Protein<REAL> proteinCenter,
//		d_Protein<REAL> proteinPartner,
//		unsigned const idx_proteinCenter,
//		unsigned const idx_proteinPartner,
//		REAL* buffer_defoX,
//		REAL* buffer_defoY,
//		REAL* buffer_defoZ,
//		REAL* buffer_trafoX,
//		REAL* buffer_trafoY,
//		REAL* buffer_trafoZ
//		)
//{
//	//using real4_t = typename TypeWrapper<REAL>::real4_t;
//	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
//	if (idx < numAtoms*numDOFs) {
//		unsigned DOFidx = idx / numAtoms;
//		auto dof = dofs[DOFidx];
//		unsigned atomIdx = idx % numAtoms;
//		Vec3<REAL> posAtomCenter(	buffer_defoX[idx],
//									buffer_defoY[idx],
//									buffer_defoZ[idx]);
//		//Vec3<REAL> deformVec(0.0f);
//
////		deform< REAL, DOF_T>(	 dof,posAtomCenter, proteinCenter, atomIdx, idx_proteinCenter,
////				deformVec.x, deformVec.y, deformVec.z);
//
//
//		device_transform(dof, posAtomCenter,proteinCenter.pivot, proteinPartner.pivot, idx_proteinCenter, idx_proteinPartner,
//		       		 buffer_trafoX[idx],  buffer_trafoY[idx], buffer_trafoZ[idx]);
//
//		//DEBUG
////		buffer_trafoX[idx] += pivotCenter.x;
////		buffer_trafoY[idx] += pivotCenter.y;
////		buffer_trafoZ[idx] += pivotCenter.z;
//
//	}
//}





template<typename REAL>
__global__ void potForce_kernel(
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		d_Protein<REAL> const  protein,
		unsigned const numDOFs,
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

		float4 potForce{0,0,0,0};
			PotForce_device( inner, outer, protein, numDOFs, idx, buffer_trafoX[idx], buffer_trafoY[idx], buffer_trafoZ[idx], potForce );


		data_out_x[idx] = potForce.x;
		data_out_y[idx] = potForce.y;
		data_out_z[idx] = potForce.z;
		data_out_E[idx] = potForce.w;
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
		double const radius_cutoff,
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
		float4 potForce{0,0,0,0};
		d_DOFPos_device< REAL, DOF_T>( protein, dof, idx, type_protein,
				 buffer_defoX[idx],
				 buffer_defoY[idx],
				 buffer_defoZ[idx],
				 x_trafo,
				 y_trafo,
				 z_trafo
				);
		if (!(radius_cutoff < 10000)){
			PotForce_device( inner, outer, protein, numDOFs, idx, x_trafo, y_trafo, z_trafo, potForce );
		}
		buffer_trafoX[idx] = x_trafo;
		buffer_trafoY[idx] = y_trafo;
		buffer_trafoZ[idx] = z_trafo;

		data_out_x[idx] = potForce.x;
		data_out_y[idx] = potForce.y;
		data_out_z[idx] = potForce.z;
		data_out_E[idx] = potForce.w;
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
		double const radius_cutoff,
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
			radius_cutoff,
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
		double const radius_cutoff,
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
		double const radius_cutoff,
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
		double const radius_cutoff,
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
		double const radius_cutoff,
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





template<typename REAL, typename DOF_T>
 void d_deform(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<REAL> const  protein,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const idx_protein,
		REAL* posX,
		REAL* posY,
		REAL* posZ,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ)
{
	cudaVerifyKernel((
			deform_kernel<<<gridSize, blockSize, 0, stream>>> (
			protein,
			dofs,
			numDOFs,
			idx_protein,
			posX,
			posY,
			posZ,
			buffer_defoX,
			buffer_defoY,
			buffer_defoZ))
		);
}

template
 void d_deform<double, DOF_MB_Modes<double>>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<double> const  protein,
		DOF_MB_Modes<double> const * dofs,
		unsigned const numDOFs,
		unsigned const idx_protein,
		double* posX,
		double* posY,
		double* posZ,
		double* buffer_defoX,
		double* buffer_defoY,
		double* buffer_defoZ);

template
 void d_deform<float, DOF_MB_Modes<float>>(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		d_Protein<float> const  protein,
		DOF_MB_Modes<float> const * dofs,
		unsigned const numDOFs,
		unsigned const idx_protein,
		float* posX,
		float* posY,
		float* posZ,
		float* buffer_defoX,
		float* buffer_defoY,
		float* buffer_defoZ);

template<typename REAL, typename DOF_T>
void d_transform(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		DOF_T const * dofs,
		unsigned const numDOFs,
		unsigned const numAtoms,
		d_Protein<REAL> const  proteinCenter,
		d_Protein<REAL> const  proteinPartner,
		unsigned const idx_proteinCenter,
		unsigned const idx_proteinPartner,
//		Vec3<REAL> pivotCenter,
//		Vec3<REAL> pivotPartner,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ
		)
{
	cudaVerifyKernel((
			transform_kernel<<<gridSize, blockSize, 0, stream>>> (
			dofs,
			numDOFs,
			numAtoms,
			proteinCenter,
			proteinPartner,
			idx_proteinCenter,
			idx_proteinPartner,
//			pivotCenter,
//			pivotPartner,
			buffer_defoX,
			buffer_defoY,
			buffer_defoZ,
			buffer_trafoX,
			buffer_trafoY,
			buffer_trafoZ
					))
		);
}

template
void d_transform<float, DOF_MB_Modes<float>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		DOF_MB_Modes<float> const * dofs,
		unsigned const numDOFs,
		unsigned const numAtoms,

//		Vec3<float> pivotCenter,
//		Vec3<float> pivotPartner,
		d_Protein<float> proteinCenter,
		d_Protein<float> proteinParnter,
				unsigned const idx_proteinCenter,
				unsigned const idx_proteinPartner,
		float* buffer_defoX,
		float* buffer_defoY,
		float* buffer_defoZ,
		float* buffer_trafoX,
		float* buffer_trafoY,
		float* buffer_trafoZ
		);

template
void d_transform<double, DOF_MB_Modes<double>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		DOF_MB_Modes<double> const * dofs,
		unsigned const numDOFs,
		unsigned const numAtoms,
//		unsigned const idx_proteinCenter,
//		unsigned const idx_proteinPartner,
//		Vec3<double> pivotCenter,
//		Vec3<double> pivotPartner,
		d_Protein<double> proteinCenter,
		d_Protein<double> proteinParnter,
		unsigned const idx_proteinCenter,
		unsigned const idx_proteinPartner,
		double* buffer_defoX,
		double* buffer_defoY,
		double* buffer_defoZ,
		double* buffer_trafoX,
		double* buffer_trafoY,
		double* buffer_trafoZ
		);

template<typename REAL>
 void d_PotForce(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<REAL> inner,
		const d_IntrlpGrid<REAL> outer,
		d_Protein<REAL> const  protein,
		unsigned const numDOFs,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ,
		REAL* data_out_x,
		REAL* data_out_y,
		REAL* data_out_z,
		REAL* data_out_E)
{
	cudaVerifyKernel((
			potForce_kernel<<<gridSize, blockSize, 0, stream>>> (
				inner,
				outer,
				protein,
				numDOFs,
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
 void d_PotForce<float>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<float> inner,
		const d_IntrlpGrid<float> outer,
		d_Protein<float> const  protein,
		unsigned const numDOFs,
		float* buffer_trafoX,
		float* buffer_trafoY,
		float* buffer_trafoZ,
		float* data_out_x,
		float* data_out_y,
		float* data_out_z,
		float* data_out_E);
template
 void d_PotForce<double>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_IntrlpGrid<double> inner,
		const d_IntrlpGrid<double> outer,
		d_Protein<double> const  protein,
		unsigned const numDOFs,
		double* buffer_trafoX,
		double* buffer_trafoY,
		double* buffer_trafoZ,
		double* data_out_x,
		double* data_out_y,
		double* data_out_z,
		double* data_out_E);



template<typename REAL>
void d_transform_new(
		 unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		dofMB<REAL>* dofSupport,
		unsigned const numDOFs,
		unsigned const numAtoms,
		unsigned const numProtein,
		unsigned const idx_proteinCenter,
		unsigned const idx_proteinPartner,
		REAL* buffer_defoX,
		REAL* buffer_defoY,
		REAL* buffer_defoZ,
		REAL* buffer_trafoX,
		REAL* buffer_trafoY,
		REAL* buffer_trafoZ
		)
{
	cudaVerifyKernel((
			transform_kernel_new<<<gridSize, blockSize, 0, stream>>> (
					dofSupport,
						numDOFs,
						numAtoms,
						numProtein,
						idx_proteinCenter,
						idx_proteinPartner,
						buffer_defoX,
						buffer_defoY,
						buffer_defoZ,
						buffer_trafoX,
						buffer_trafoY,
						buffer_trafoZ
					))
		);
}

template
void d_transform_new<float>(
		 unsigned blockSize,
				unsigned gridSize,
				const cudaStream_t &stream,
		dofMB<float>* dofSupport,
			unsigned const numDOFs,
			unsigned const numAtoms,
			unsigned const numProtein,
			unsigned const idx_proteinCenter,
			unsigned const idx_proteinPartner,
			float* buffer_defoX,
			float* buffer_defoY,
			float* buffer_defoZ,
			float* buffer_trafoX,
			float* buffer_trafoY,
			float* buffer_trafoZ
			);

template
void d_transform_new<double>(
		 unsigned blockSize,
				unsigned gridSize,
				const cudaStream_t &stream,
		dofMB<double>* dofSupport,
			unsigned const numDOFs,
			unsigned const numAtoms,
			unsigned const numProtein,
			unsigned const idx_proteinCenter,
			unsigned const idx_proteinPartner,
			double* buffer_defoX,
			double* buffer_defoY,
			double* buffer_defoZ,
			double* buffer_trafoX,
			double* buffer_trafoY,
			double* buffer_trafoZ
			);

}
