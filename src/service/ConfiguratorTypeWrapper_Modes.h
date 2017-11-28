/*
 * ConfiguratorTypeWrapper_Modes.h
 *
 *  Created on: Aug 15, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATORTYPEWRAPPER_MODES_H_
#define SRC_SERVICE_CONFIGURATORTYPEWRAPPER_MODES_H_

#include <type_traits>

#include "Configurator_6D_Modes.h"
#include "Types_6D_Modes.h"

namespace as {

template<typename GenericTypes>
struct ConfiguratorTypeWrapper_Modes{
	ConfiguratorTypeWrapper_Modes() = delete;

	using configurator_t_Modes = typename std::conditional<std::is_same<GenericTypes, Types_6D_Modes<float>>::value,
			Configurator_6D_Modes<float>, Configurator_6D_Modes<double>>::type;

};

} // namespace



#endif /* SRC_SERVICE_CONFIGURATORTYPEWRAPPER_H_ */
