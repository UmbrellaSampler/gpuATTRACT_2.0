/*
 * transform_modes.h
 *
 *  Created on: Dec 5, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TRANSFORM_MODES_H_
#define SRC_MODEL_TRANSFORM_MODES_H_


#ifdef CUDA

template<typename REAL>
void d_DOFPos(
	unsigned blockSize,
	unsigned gridSize,
	const cudaStream_t &stream,
	REAL const* xRec,
	REAL const* yRec,
	REAL const* zRec,
	REAL const* xLig,
	REAL const* yLig,
	REAL const* zLig,
	REAL const* xModesRec,
	REAL const* yModesRec,
	REAL const* zModesRec,
	REAL const* xModesLig,
	REAL const* yModesLig,
	REAL const* zModesLig,
	DOF_6D_Modes<REAL>* dofs,
	unsigned numAtomsRec,
	unsigned numAtomsLig,
	unsigned numModesRec,
	unsigned numModesLig,
	unsigned numDOFsLig,
	REAL* xRecDefo,
	REAL* yRecDefo,
	REAL* zRecDefo,
	REAL* xRecTrafo,
	REAL* yRecTrafo,
	REAL* zRecTrafo,
	REAL* xLigTrafo,
	REAL* yLigTrafo,
	REAL* zLigTrafo);

#endif


#endif /* TRANSFORM_MODES_H_ */
