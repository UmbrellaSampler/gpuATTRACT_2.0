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

namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_6D_Impl;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_6D_Impl<REAL> const&);

template<typename REAL>
struct Result_6D_Impl;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_6D_Impl<REAL> const& args);

#endif

template<typename REAL>
struct DOF_6D_Impl {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	vec3_t pos;
	vec3_t ang;
};

template<typename REAL>
struct Result_6D_Impl {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	real_t E;
	vec3_t pos;
	vec3_t ang;
};

template<typename REAL>
struct Types_6D {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;

	using DOF = DOF_6D_Impl<real_t>;

	struct Common {
		id_t gridId;
		id_t ligId;
		id_t recId;
		id_t tableId;
		id_t paramsId;
	};

	using Result = Result_6D_Impl<REAL>;
};

}  // namespace as



#endif /* SRC_TYPES_6D_H_ */
