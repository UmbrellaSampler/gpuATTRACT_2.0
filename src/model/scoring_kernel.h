/*
 * scoring_kernel.h
 *
 *  Created on: Jun 12, 2018
 *      Author: glenn
 */

#ifndef SCORING_KERNEL_H_
#define SCORING_KERNEL_H_
#include "transform_modes.h"
#include "Protein.h"
#include "IntrplGrid.h"
#include "nativeTypesFunctions.h"
#include "nativeTypesMath.h"
#include "VoxelOctet.h"
#include "trilinIntrpl.h"
#include <cmath>
#include "Types_MB_Modes.h"

#ifdef CUDA
#include "nativeTypesWrapper.h"
#include "DeviceIntrplGrid.h"
#include "DeviceProtein.h"
#include "macros.h"
#endif
namespace as{

template <typename REAL>
void getNewDof(unsigned numDOFs, unsigned& numProtein, DOF_MB_Modes<REAL> const * dofs, std::vector<Vec3<REAL>>& pivots, dofMB<REAL> * dofs_out){

//	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
//
//	const unsigned idx_protein =  idx % numProtein;
//	const unsigned idx_dof = idx / numProtein;
//	if (idx < numDOFs* numProtein) {
	for( int id_dof = 0; id_dof< numDOFs ; ++id_dof){
		for( int id_center = 0; id_center<numProtein ; ++id_center){
			for ( unsigned id_partner = 0; id_partner < numProtein; ++ id_partner){
				if ( id_partner != id_center){
					auto const dof = dofs[id_dof];
					auto const proteinActive = dof.protein[id_center];

					RotMat<REAL>  rotMatCenter = euler2rotmat(proteinActive.ang.x,
															  proteinActive.ang.y,
															  proteinActive.ang.z);
					RotMat<REAL>  rotMatPartner = euler2rotmat(dof.protein[id_partner].ang.x,
															   dof.protein[id_partner].ang.y,
															   dof.protein[id_partner].ang.z).getInv();
		//			REAL rotMat[9] = {0};//= rotMatPartner * rotMatCenter;
		//			for(unsigned i = 0; i < 3; ++i) {
		//				for(unsigned j = 0; j < 3; ++j) {
		//					for(unsigned k = 0; k < 3; ++k) {
		//						rotMat[i*3 + j] += rotMatCenter[i*3 + k]*rotMatPartner[k*3 + j];
		//					}
		//				}
		//			}
					RotMat<REAL> rotMat =  rotMatPartner *rotMatCenter;

		//			dofMB<REAL> dofMB;
		//			dofMB.rotMat = rotMat;
		//			dofMB.translation = rotMatPartner * (proteinActive.pos - dof.protein[id_partner].pos + pivots[idx_protein] -  pivots[id_partner]);
					//rotMatPartner * (proteinActive.pos - dof.protein[id_partner].pos + pivots[idx_protein] -  pivots[id_partner]);

//					std::cout << "pivot "<<id_center << id_partner << pivots[id_center] << pivots[id_partner]<< std::endl;
//					std::cout << "pivot "<<id_center << id_partner << proteinActive.pos << dof.protein[id_partner].pos<< std::endl;
					dofs_out[numDOFs * numProtein * id_center + numProtein * id_dof + id_partner].translation = rotMatPartner * (proteinActive.pos - dof.protein[id_partner].pos + pivots[id_center] -  pivots[id_partner]);
					dofs_out[numDOFs * numProtein * id_center + numProtein * id_dof + id_partner].rotMat = rotMat;
					//std::cout <<"pos"<< proteinActive.pos << dof.protein[id_partner].pos << std::endl;
//					std::cout << id_center << id_partner << dofs_out[numDOFs * numProtein * id_center + numProtein * id_dof + id_partner].translation << std::endl;
//					std::cout << id_center << id_partner << dofs_out[numDOFs * numProtein * id_center + numProtein * id_dof + id_partner].rotMat << std::endl;
				}
			}
		}
	}
}

template
void getNewDof<float>(unsigned numDOFs, unsigned& numProtein, DOF_MB_Modes<float> const* dofs, std::vector<Vec3<float>>& pivots, dofMB<float> * dofs_out);



template
void getNewDof<double>(unsigned numDOFs, unsigned& numProtein, DOF_MB_Modes<double> const * dofs, std::vector<Vec3<double>>& pivots, dofMB<double> * dofs_out);


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
		REAL* data_out_E);

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
		REAL* buffer_defoZ);



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
		REAL* data_out_E);

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
		);
}
#endif /* SCORING_KERNEL_H_ */
