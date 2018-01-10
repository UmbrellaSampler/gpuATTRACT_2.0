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
#include "Configurator_6D_Modes.h"
#include "Configurator_6D_MB_Modes.h"
#include "Types_6D.h"
#include "Types_6D_Modes.h"
#include "Types_6D_MB_Modes.h"

namespace as {

template<typename GenericTypes>
struct ConfiguratorTypeWrapper {
	ConfiguratorTypeWrapper() = delete;

	using configurator_t =
			typename std::conditional<std::is_same<GenericTypes, Types_6D<float>>::value,
			Configurator_6D<float>,
			typename std::conditional<std::is_same<GenericTypes, Types_6D<double>>::value,
			Configurator_6D<double>,
			typename std::conditional<std::is_same<GenericTypes, Types_6D_Modes<float>>::value,
			Configurator_6D_Modes<float>,
			typename std::conditional<std::is_same<GenericTypes, Types_6D_Modes<double>>::value,
			Configurator_6D_Modes<double>,
			typename std::conditional<std::is_same<GenericTypes, Types_6D_MB_Modes<float>>::value,
			Configurator_6D_MB_Modes<float>,
			Configurator_6D_MB_Modes<double>>
			::type>
			::type>
			::type>
			::type>
			::type;

};

} // namespace



#endif /* SRC_SERVICE_CONFIGURATORTYPEWRAPPER_H_ */
