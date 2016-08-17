/*
 * DOF.h
 *
 *  Created on: Aug 14, 2016
 *      Author: uwe
 */

#ifndef SRC_DOF_6D_H_
#define SRC_DOF_6D_H_

#include "nativeTypesWrapper.h"
#include "Vec3.h"

namespace as {

template<typename REAL>
class DOF_6D {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
public:
	vec3_t pos;
	vec3_t ang;
};

} // namespace



#endif /* SRC_DOF_H_ */
