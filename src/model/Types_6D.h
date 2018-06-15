/*
5 * Types_2B_6D.h
 *
 *  Created on: Aug 9, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_6D_H_
#define SRC_TYPES_6D_H_

#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "GenericTypes.h"

namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_6D;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_6D<REAL> const&);

template<typename REAL>
struct Result_6D;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_6D<REAL> const& args);

#endif

template<typename REAL>
struct DOF_6D {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	vec3_t pos;
	vec3_t ang;
};


struct Common {
	using real_t = double;
	using vec3_t = Vec3<real_t>;
	id_t gridId;
	id_t ligId;
	id_t recId;
	id_t tableId;
	id_t paramsId;
	vec3_t pivotRec;
	vec3_t pivotLig;
	real_t radius_cutoff;
};

template<typename REAL>
struct Result_6D {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	real_t E;
	vec3_t pos;
	vec3_t ang;
	REAL get_Energy(){ return E;}
};

template<typename REAL>
using Types_6D = GenericTypes<DOF_6D<REAL>, Common, Result_6D<REAL>>;

}  // namespace as



#endif /* SRC_TYPES_6D_H_ */
