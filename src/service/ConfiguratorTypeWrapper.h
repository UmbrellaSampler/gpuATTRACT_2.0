/*
 * ConfiguratorTypeWrapper.h
 *
 *  Created on: Aug 15, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATORTYPEWRAPPER_H_
#define SRC_SERVICE_CONFIGURATORTYPEWRAPPER_H_

#include <type_traits>

#include "Configurator_6D.h"
#include "Types_6D.h"

namespace as {

//template<typename REAL, template <typename REAL> class GenericTypes>
template<typename GenericTypes>
struct ConfiguratorTypeWrapper {
//	static_assert(std::is_arithmetic<REAL>::value, "Only arithmetic types supported");
	ConfiguratorTypeWrapper() = delete;

	using configurator_t = typename std::conditional<std::is_same<GenericTypes, Types_6D<float>>::value,
			Configurator_6D<float>, Configurator_6D<double>>::type;


};

} // namespace



#endif /* SRC_SERVICE_CONFIGURATORTYPEWRAPPER_H_ */
