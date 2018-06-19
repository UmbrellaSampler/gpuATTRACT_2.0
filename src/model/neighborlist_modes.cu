/*
 * neighborlist.cu
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#include "neighborlist_modes.h"

#include "nativeTypesWrapper.h"
#include "DeviceNLGrid.h"
#include "DeviceProtein.h"
#include "DeviceParamTable.h"
#include "SimParam.h"
#include "forcefield.h"
#include "macros.h"
#include "Types_6D_Modes.h"
#include "transform_modes.h"

namespace as {

/*
 * In comparison to the NL force calculation without modes, this methode takes the deformed coordinates of the receptor instead of its original position,
 * which does not contain mode deformation.
 * Note that if the forces that are acting on the receptor are calculated, they have to be rotated back into the system of the receptor.
 */
template<typename REAL, typename DOF_T>
__global__ void d_NLPotForce(
		const d_NLGrid<REAL> grid,
		const d_Protein<REAL> rec,
		const d_Protein<REAL> lig,
		const d_ParamTable<REAL> table,
		const DOF_T* dofs,
		const double radius_cutoff,
		const SimParam<REAL> simParam,
		const unsigned numDOFs,
		const REAL* RecPosX,
		const REAL* RecPosY,
		const REAL* RecPosZ,
		const REAL* LigPosX,
		const REAL* LigPosY,
		const REAL* LigPosZ,
		REAL* outLig_fx,
		REAL* outLig_fy,
		REAL* outLig_fz,
		REAL* outLigand_E)
{
	using real3_t = typename TypeWrapper<REAL>::real3_t;
	const unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned LigNumEl = lig.numAtoms;
	if (i < LigNumEl*numDOFs) {
		const unsigned LigAttrIdx = i % LigNumEl;
		unsigned recBase = 0;


		const unsigned atomTypeLig = lig.type[LigAttrIdx];

		if (atomTypeLig != 0) {


			const REAL posLigX = LigPosX[i];
			const REAL posLigY = LigPosY[i];
			const REAL posLigZ = LigPosZ[i];

			/* test if particle is out of bounds and perform data fetch and neigbourlist calculations */
			if (!(     (posLigX < grid.minDim.x || posLigX > grid.maxDim.x)
					|| (posLigY < grid.minDim.y || posLigY > grid.maxDim.y)
					|| (posLigZ < grid.minDim.z || posLigZ > grid.maxDim.z) ))
			{

				const uint2 nDesc = tex3D<uint2>(grid.tex,
						(posLigX - grid.minDim.x) * grid.dVox_inv + 0.5,
						(posLigY - grid.minDim.y) * grid.dVox_inv + 0.5,
						(posLigZ - grid.minDim.z) * grid.dVox_inv + 0.5);
				/* numEl = x; idx = y */


				real3_t fAcc = {0,0,0};
				REAL eAcc = 0;
				for (unsigned j = 0; j < nDesc.x; ++j) {
					const unsigned nIdx = grid.neighborList[nDesc.y + j];

					REAL xRec, yRec, zRec;
					if( std::is_same<DOF_T, DOF_6D_Modes<REAL>>::value ){
						recBase = rec.numAtoms * (int) (i / LigNumEl);
						xRec = RecPosX[nIdx + recBase];
						yRec = RecPosY[nIdx + recBase];
						zRec = RecPosZ[nIdx + recBase];
					}
					else{
						xRec = rec.xPos[nIdx];
						yRec = rec.yPos[nIdx];
						zRec = rec.zPos[nIdx];
					}


					REAL dx = posLigX - xRec;
					REAL dy = posLigY - yRec;
					REAL dz = posLigZ - zRec;

					const REAL dr2 = dx * dx + dy * dy + dz * dz;
					const REAL dPlateau2 = grid.dPlateau2;
					if ((dr2) > dPlateau2) {
						continue;
					}

					constexpr REAL one = static_cast<REAL>(1.0);
					const REAL dr2_inv = one/dr2; // inverse of dr2

					// Scale distances
					dx *= dr2_inv;
					dy *= dr2_inv;
					dz *= dr2_inv;

					real3_t fVdW;
					REAL eVdW;

					const size_t atomTypeRec = rec.type[nIdx];


					// calculate energy and potential/energy of LJ/VdW potential

					auto const params = table.getParams(atomTypeRec-1, atomTypeLig-1);
					LJPotForce(dr2, dr2_inv, dx, dy, dz,
							params,
							one, table.shape,
							fVdW.x, fVdW.y, fVdW.z, eVdW);

					fAcc.x  += fVdW.x;
					fAcc.y  += fVdW.y;
					fAcc.z  += fVdW.z;
					eAcc += eVdW;

					const REAL chargeLig = lig.charge[LigAttrIdx];
					const REAL chargeRec = rec.charge[nIdx];
					const REAL chargeLigRec = chargeLig * chargeRec * simParam.ffelec;

					const bool calc_elec = abs(chargeLigRec) > 0.001; // evaluate electric potential

					REAL dPlateau2_inv = 1/grid.dPlateau2;
					const REAL ratio = sqrt(dr2*dPlateau2_inv);
					REAL rdx = ratio*dx;
					REAL rdy = ratio*dy;
					REAL rdz = ratio*dz;

					LJPotForce(dPlateau2, dPlateau2_inv, rdx, rdy, rdz,
						params,
						one, table.shape,
						fVdW.x, fVdW.y, fVdW.z, eVdW);
					fAcc.x  -= fVdW.x;
					fAcc.y  -= fVdW.y;
					fAcc.z  -= fVdW.z;
					eAcc -= eVdW;


					if (calc_elec) {
						REAL eEl;
						real3_t fEl;

						// calculate energy and potential/energy of charge potential

						if (false) {
							printf("%u %f %f %f %u\n" ,
									i, posLigX, posLigY, posLigZ, atomTypeLig);
						}

						ChargePotForce(dr2_inv, dx, dy, dz,
								chargeLigRec,
								one, simParam.dielec,
								fEl.x, fEl.y, fEl.z, eEl);

						fAcc.x += fEl.x;
						fAcc.y += fEl.y;
						fAcc.z += fEl.z;
						eAcc += eEl;

						ChargePotForce(dPlateau2_inv, rdx, rdy, rdz,
								chargeLigRec,
								one, simParam.dielec,
								fEl.x, fEl.y, fEl.z, eEl);
						fAcc.x -= fEl.x;
						fAcc.y -= fEl.y;
						fAcc.z -= fEl.z;
						eAcc -= eEl;

					}
				}


				/* store results back to global memory */
				if (nDesc.x > 0) {
					fAcc.x = fAcc.x + outLig_fx[i];
					fAcc.y = fAcc.y + outLig_fy[i];
					fAcc.z = fAcc.z + outLig_fz[i];
					if(lig.isOrigin ){
						device_rotateForces<REAL, DOF_T>( fAcc.x, fAcc.y, fAcc.z, dofs[i/LigNumEl] );
					}
					outLig_fx[i] = fAcc.x;
					outLig_fy[i] = fAcc.y;
					outLig_fz[i] = fAcc.z;

					outLigand_E[i] += eAcc;
				}
			}
		} // if (atomtype != 0)
	}
}




template<typename REAL, typename DOF_T>
void d_NLPotForce(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_NLGrid<REAL>& grid,
		const d_Protein<REAL>& rec,
		const d_Protein<REAL>& lig,
		const d_ParamTable<REAL>& table,
		const DOF_T* dofs,
		const double radius_cutoff,
		const SimParam<REAL>& simParam,
		const unsigned& numDOFs,
		const REAL* RecPosX,
		const REAL* RecPosY,
		const REAL* RecPosZ,
		const REAL* LigPosX,
		const REAL* LigPosY,
		const REAL* LigPosZ,
		REAL* outLig_fx,
		REAL* outLig_fy,
		REAL* outLig_fz,
		REAL* outLigand_E)
{
	cudaVerifyKernel((
			d_NLPotForce<REAL,DOF_T><<<gridSize, blockSize, 0, stream>>> (
				grid,
				rec,
				lig,
				table,
				dofs,
				radius_cutoff,
				simParam,
				numDOFs,
				RecPosX,
				RecPosY,
				RecPosZ,
				LigPosX,
				LigPosY,
				LigPosZ,
				outLig_fx,
				outLig_fy,
				outLig_fz,
				outLigand_E
			)
		));
}

template
void d_NLPotForce<float, DOF_6D_Modes<float>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_NLGrid<float>& grid,
		const d_Protein<float>& rec,
		const d_Protein<float>& lig,
		const d_ParamTable<float>& table,
		const DOF_6D_Modes<float>* dofs,
		const double radius_cutoff,
		const SimParam<float>& simParam,
		const unsigned& numDOFs,
		const float* RecPosX,
		const float* RecPosY,
		const float* RecPosZ,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		float* outLig_fx,
		float* outLig_fy,
		float* outLig_fz,
		float* outLigand_E
		);

template
void d_NLPotForce<double, DOF_6D_Modes<double>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_NLGrid<double>& grid,
		const d_Protein<double>& rec,
		const d_Protein<double>& lig,
		const d_ParamTable<double>& table,
		const DOF_6D_Modes<double>* dofs,
		const double radius_cutoff,
		const SimParam<double>& simParam,
		const unsigned& numDOFs,
		const double* RecPosX,
		const double* RecPosY,
		const double* RecPosZ,
		const double* LigPosX,
		const double* LigPosY,
		const double* LigPosZ,
		double* outLig_fx,
		double* outLig_fy,
		double* outLig_fz,
		double* outLigand_E
		);

template
void d_NLPotForce<float, DOF_6D<float>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_NLGrid<float>& grid,
		const d_Protein<float>& rec,
		const d_Protein<float>& lig,
		const d_ParamTable<float>& table,
		const DOF_6D<float>* dofs,
		const double radius_cutoff,
		const SimParam<float>& simParam,
		const unsigned& numDOFs,
		const float* RecPosX,
		const float* RecPosY,
		const float* RecPosZ,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		float* outLig_fx,
		float* outLig_fy,
		float* outLig_fz,
		float* outLigand_E
		);

template
void d_NLPotForce<double, DOF_6D<double>>(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_NLGrid<double>& grid,
		const d_Protein<double>& rec,
		const d_Protein<double>& lig,
		const d_ParamTable<double>& table,
		const DOF_6D<double>* dofs,
		const double radius_cutoff,
		const SimParam<double>& simParam,
		const unsigned& numDOFs,
		const double* RecPosX,
		const double* RecPosY,
		const double* RecPosZ,
		const double* LigPosX,
		const double* LigPosY,
		const double* LigPosZ,
		double* outLig_fx,
		double* outLig_fy,
		double* outLig_fz,
		double* outLigand_E
		);


}  // namespace as

