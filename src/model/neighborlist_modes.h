/*
 * neighborlist_modes.h
 *
 *  Created on: Dec 4, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_NEIGHBORLIST_MODES_H_
#define SRC_MODEL_NEIGHBORLIST_MODES_H_


#include "neighborlist.h"


namespace as {

/**
 * this is an overloaded Version of the original  NLPotForce function which directly calculates the counteracting forces acting on the receptor.
 */
template<typename REAL>
void NLPotForce(
		NLGrid<REAL> const* grid,
		Protein<REAL> const* rec, //TODO: use receptor positions instead of receptor
		Protein<REAL> const* lig,
		SimParam<REAL> const* simParam,
		ParamTable<REAL> const* table,
		REAL const* RecPosX,
		REAL const* RecPosY,
		REAL const* RecPosZ,
		REAL const* LigPosX,
		REAL const* LigPosY,
		REAL const* LigPosZ,
		REAL* outLig_fx,
		REAL* outLig_fy,
		REAL* outLig_fz,
		REAL* outLig_E,
		REAL* outRec_fx,// receptor forces have to be set to zero in advance
		REAL* outRec_fy,
		REAL* outRec_fz
		)
{

	using real3_t = typename TypeWrapper<REAL>::real3_t;
	using real4_t = typename TypeWrapper<REAL>::real4_t;

	const unsigned numAtomsLig = lig->numAtoms();
	for (unsigned i = 0; i < numAtomsLig; ++i) {
		const unsigned atomTypeLig = lig->type()[i];

		if (atomTypeLig == 0)
			continue;

		const REAL posLigX = LigPosX[i];
		const REAL posLigY = LigPosY[i];
		const REAL posLigZ = LigPosZ[i];


		/* test if particle is out of bounds and perform data fetch and neighbourlist calculations */
		if (!(grid->outOfBounds(posLigX, posLigY, posLigZ))) {

			int idxX, idxY, idxZ;
			grid->getIndex(posLigX, posLigY, posLigZ, idxX, idxY, idxZ);
			const NeighbourDesc &nDesc = grid->getNeighbourDesc(idxX, idxY,
					idxZ);

			real3_t fRec = {0.0,0.0,0.0};
			real3_t fAcc = {0.0,0.0,0.0};
			REAL eAcc = 0.0;




			for (unsigned j = 0; j < nDesc.numEl; ++j) {
				const unsigned nIdx = grid->getNeighbor(nDesc.idx + j);

				REAL dx = posLigX - RecPosX[nIdx];
				REAL dy = posLigY - RecPosY[nIdx];
				REAL dz = posLigZ - RecPosZ[nIdx];

				fRec.x = 0.0;
				fRec.y = 0.0;
				fRec.z = 0.0;



				const REAL dr2 = dx * dx + dy * dy + dz * dz;

				if (grid->outOfPlateau(dr2)) {
					continue;
				}

				const REAL dr2_inv = 1.0/dr2; // inverse of dr2

				// Scale distances
				dx *= dr2_inv;
				dy *= dr2_inv;
				dz *= dr2_inv;

				const size_t atomTypeRec = rec->type()[nIdx];

				// calculate energy and potential/energy of LJ/VdW potential
				assert(atomTypeRec > 0);
				assert(atomTypeRec < 99);
				assert(atomTypeLig > 0);
				assert(atomTypeLig < 99);

				real3_t fVdW;
				REAL eVdW;

				REAL swi = 1;
				REAL swiPlateau = 1;
				LJPotForce(dr2, dr2_inv, dx, dy, dz,
						table->getParams(atomTypeRec-1, atomTypeLig-1),
						swi, table->potShape(),
						fVdW.x, fVdW.y, fVdW.z, eVdW);

				fAcc.x  += fVdW.x;
				fAcc.y  += fVdW.y;
				fAcc.z  += fVdW.z;
				eAcc += eVdW;

				fRec.x -= fVdW.x;
				fRec.y -= fVdW.y;
				fRec.z -= fVdW.z;


				const REAL chargeLig = lig->charge()[i];
				const REAL chargeRec = rec->charge()[nIdx];

				const REAL chargeLigRec = chargeLig * chargeRec * simParam->ffelec;
				const bool calc_elec = abs(chargeLigRec) > 0.001; // evaluate electric potential

				REAL rdx, rdy, rdz;
				REAL ratio = grid->getRatio(dr2);
				rdx = ratio*dx;
				rdy = ratio*dy;
				rdz = ratio*dz;

				LJPotForce(grid->dPlateau2(), grid->dPlateau2_inv(), rdx, rdy, rdz,
					table->getParams(atomTypeRec-1, atomTypeLig-1),
					swiPlateau, table->potShape(),
					fVdW.x, fVdW.y, fVdW.z, eVdW);

				fAcc.x  -= fVdW.x;
				fAcc.y  -= fVdW.y;
				fAcc.z  -= fVdW.z;
				eAcc -= eVdW;

				fRec.x += fVdW.x;
				fRec.y += fVdW.y;
				fRec.z += fVdW.z;



				if (calc_elec) {
					REAL eEl;
					real4_t fEl;

					// calculate energy and potential/energy of charge potential
					ChargePotForce(dr2_inv, dx, dy, dz,
							chargeLigRec,
							swi, simParam->dielec,
							fEl.x, fEl.y, fEl.z, eEl);

					fAcc.x += fEl.x;
					fAcc.y += fEl.y;
					fAcc.z += fEl.z;
					eAcc += eEl;

					fRec.x -= fEl.x;;
					fRec.y -= fEl.y;
					fRec.z -= fEl.z;


					ChargePotForce(grid->dPlateau2_inv(), rdx, rdy, rdz,
							chargeLigRec,
							swiPlateau, simParam->dielec,
							fEl.x, fEl.y, fEl.z, eEl);
					fAcc.x -= fEl.x;
					fAcc.y -= fEl.y;
					fAcc.z -= fEl.z;
					eAcc -= eEl;

					fRec.x += fEl.x;;
					fRec.y += fEl.y;
					fRec.z += fEl.z;



				}

				outRec_fx[nIdx] += fRec.x;
				outRec_fy[nIdx] += fRec.y;
				outRec_fz[nIdx] += fRec.z;
			} // for j

			outLig_fx[i] += fAcc.x;
			outLig_fy[i] += fAcc.y;
			outLig_fz[i] += fAcc.z;
			outLig_E[i]  += eAcc;



		} // for i
	}
}

#ifdef CUDA

template<typename REAL>
void d_NLPotForce(
		unsigned blockSize,
		unsigned gridSize,
		const cudaStream_t &stream,
		const d_NLGrid<REAL>& grid,
		const d_Protein<REAL>& rec,
		const d_Protein<REAL>& lig,
		const d_ParamTable<REAL>& table,
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
		REAL* outLigand_E);

#endif

}//end Namespace as

#endif /* NEIGHBORLIST_MODES_H_ */
