/*
 * scATTRACT.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#ifndef SRC_SCATTRACT_H_
#define SRC_SCATTRACT_H_

#include <memory>
#include "App.h"
#include "ConfiguratorTypeWrapper.h"

namespace as {

template<typename GenericTypes>
class scATTRACT : public App {

	using input_t = typename GenericTypes::input_t;
	using common_t = typename GenericTypes::common_t;
	using result_t = typename GenericTypes::result_t;

	using configurator_t = typename ConfiguratorTypeWrapper<GenericTypes>::configurator_t;
	using real_t = typename configurator_t::real_t;

	using service_t = Service<GenericTypes>;
public:
	scATTRACT();
	virtual ~scATTRACT() {}

	void init(CmdArgs const& args) override;
	void finalize() override;
	void run() override;

private:

	std::shared_ptr<configurator_t> _config;

};

}  // namespace as



#endif /* SRC_SCATTRACT_H_ */
