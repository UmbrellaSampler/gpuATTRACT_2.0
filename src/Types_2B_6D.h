/*
 * Types_2B_6D.h
 *
 *  Created on: Aug 9, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_2B_6D_H_
#define SRC_TYPES_2B_6D_H_

#include "nativeTypesWrapper.h"
#include "Vec3.h"

namespace as {

template<typename REAL>
struct Types_6D {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;

	struct DOF {
		vec3_t pos;
		vec3_t ang;
	};

	struct Common {
		id_t gridId;
		id_t ligId;
		id_t recId;
		id_t tableId;
		id_t paramsId;
	};

	struct Result {
		real_t E;
		vec3_t pos;
		vec3_t ang;
	};
};

}  // namespace as



#endif /* SRC_TYPES_2B_6D_H_ */
