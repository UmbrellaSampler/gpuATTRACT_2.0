/*
 * nativeTypesWrapper.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_NATIVETYPESWRAPPER_H_
#define SRC_NATIVETYPESWRAPPER_H_

#include <type_traits>
#include "nativeTypes.h"

template<typename REAL>
struct TypeWrapper {
	static_assert(std::is_arithmetic<REAL>::value, "Only arithmetic types supported");
	TypeWrapper() = delete;

	using real_t = REAL;
	using real2_t = typename std::conditional<std::is_same<real_t, float>::value, float2, double2>::type;
	using real3_t = typename std::conditional<std::is_same<real_t, float>::value, float3, double3>::type;
	using real4_t = typename std::conditional<std::is_same<real_t, float>::value, float4, double4>::type;

};



#endif /* SRC_NATIVETYPESWRAPPER_H_ */
