/*
 * Configurator_MB_Modes.h
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATOR_MB_MODES_H_
#define SRC_SERVICE_CONFIGURATOR_MB_MODES_H_

#include <memory>
#include "CmdArgs.h"
#include "ConfiguratorBase.h"
#include "Types_MB_Modes.h"

namespace as {

template<typename REAL>
class Configurator_MB_Modes : public ConfiguratorBase<Types_MB_Modes<REAL>> {
	using genericTypes_t = Types_MB_Modes<REAL>;
	using typename ConfiguratorBase<genericTypes_t>::service_t;
	using typename ConfiguratorBase<genericTypes_t>::input_t;
	using typename ConfiguratorBase<genericTypes_t>::common_t;

public:
	using typename ConfiguratorBase<genericTypes_t>::real_t;
	using typename ConfiguratorBase<genericTypes_t>::server_t;

	void init(CmdArgs const& args) override;
	void finalize() noexcept override;

	virtual ~Configurator_MB_Modes() {}

};

}  // namespace as




#endif /* SRC_SERVICE_CONFIGURATOR_MB_MODES_H_ */
