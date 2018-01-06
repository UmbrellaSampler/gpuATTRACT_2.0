/*
 * Types_6D_modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TYPES_6D_MB_MODES_H_
#define SRC_MODEL_TYPES_6D_MB_MODES_H_

#include "Types_6D_Config.h"
#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "GenericTypes.h"
#include "Types_6D_Modes.h"

namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_6D_MB_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_6D_Modes<REAL> const& args);

template<typename REAL>
struct Result_6D_MB_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<REAL> const& args);

#endif

template<typename REAL>
struct DOF_6D_MB_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	DOF_6D<real_t> _6D[LIGANDS_MAX_NUMBER];
	real_t modesRec[MODES_MAX_NUMBER];
	real_t modesLig[LIGANDS_MAX_NUMBER][MODES_MAX_NUMBER];

};

struct Common_MB_Modes {
	id_t gridIdRec;
	id_t gridIdLig[LIGANDS_MAX_NUMBER];
	id_t ligId[LIGANDS_MAX_NUMBER];
	id_t recId;
	id_t tableId;
	id_t paramsId;
	static unsigned int numModesRec;
	static unsigned int numModesLig[LIGANDS_MAX_NUMBER];
};

template<typename REAL>
struct Result_6D_MB_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	real_t E;
	DOF_6D<real_t> _6D[LIGANDS_MAX_NUMBER];
	real_t modesRec[MODES_MAX_NUMBER];
	real_t modesLig[LIGANDS_MAX_NUMBER][MODES_MAX_NUMBER];
};

template<typename REAL>
using Types_6D_MB_Modes = GenericTypes<DOF_6D_MB_Modes<REAL>, Common_MB_Modes, Result_6D_MB_Modes<REAL>>;

}  // namespace as




#endif /* TYPES_6D_MB_MODES_H_ */
