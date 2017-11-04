/*
 * Types_6D_modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TYPES_6D_MODES_H_
#define SRC_MODEL_TYPES_6D_MODES_H_

#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "GenericTypes.h"

namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_6D_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_6D_Modes<REAL> const&);

template<typename REAL>
struct Result_6D_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<REAL> const& args);

#endif

template<typename REAL>
struct DOF_6D_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	vec3_t pos;
	vec3_t ang;
	real_t modes[10];
};

struct Common {
	id_t gridId;
	id_t ligId;
	id_t recId;
	id_t tableId;
	id_t paramsId;
};

template<typename REAL>
struct Result_6D_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	real_t E;
	vec3_t pos;
	vec3_t ang;
	real_t modes[10];
};

template<typename REAL>
using Types_6D_Modes = GenericTypes<DOF_6D_Modes<REAL>, Common, Result_6D_Modes<REAL>>;

}  // namespace as




#endif /* TYPES_6D_MODES_H_ */
