/*
 * emATTRACT.h
 *
 *  Created on: Jul 11, 2017
 *      Author: uwe
 */

#ifndef SRC_APPS_EM_EMATTRACT_H_
#define SRC_APPS_EM_EMATTRACT_H_

#include <memory>
#include "App.h"
#include "ConfiguratorTypeWrapper.h"

namespace as {

template<typename GenericTypes>
class Service;

template<typename GenericTypes>
class emATTRACT : public App {
	using dof_t = typename GenericTypes::input_t;
	using common_t = typename GenericTypes::common_t;
	using result_t = typename GenericTypes::result_t;

	using configurator_t = typename ConfiguratorTypeWrapper<GenericTypes>::configurator_t;
	using real_t = typename configurator_t::real_t;

	using service_t = Service<GenericTypes>;


	using server_t = typename configurator_t::server_t;
public:
	emATTRACT();
	virtual ~emATTRACT() {}

	void init(CmdArgs const& args) override;
	void finalize() override;
	void run() override;

private:

	std::shared_ptr<configurator_t> _config;
	std::string _solverName;

};

}  // namespace as



#endif /* SRC_APPS_EM_EMATTRACT_H_ */
