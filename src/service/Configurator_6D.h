/*
 * Init_6D.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#ifndef SRC_CONFIGURATOR_6D_H_
#define SRC_CONFIGURATOR_6D_H_

#include <memory>
#include "CmdArgs.h"
#include "ConfiguratorBase.h"
#include "Types_6D.h"
#include "Service.h"

namespace as {

template<typename REAL>
class Configurator_6D : public ConfiguratorBase<Types_6D<REAL>> {
	using genericTypes_t = Types_6D<REAL>;
	using typename ConfiguratorBase<genericTypes_t>::service_t;
	using typename ConfiguratorBase<genericTypes_t>::input_t;
	using typename ConfiguratorBase<genericTypes_t>::common_t;

public:
	using typename ConfiguratorBase<genericTypes_t>::real_t;
	using typename ConfiguratorBase<genericTypes_t>::server_t;

	void init(CmdArgs const& args) override;
	void finalize() noexcept override;

	virtual ~Configurator_6D() {}

};

}  // namespace as



#endif /* SRC_INIT_6D_H_ */
