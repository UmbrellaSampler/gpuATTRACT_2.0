/*
 * scATTRACTModes.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#ifndef SRC_SCATTRACTMODES_H_
#define SRC_SCATTRACTMODES_H_

#include <memory>
#include "App.h"
#include "ConfiguratorTypeWrapper_Modes.h"

namespace as {

template<typename GenericTypes>
class scATTRACTModes : public App {
//	using real_t = typename SERVICE::real_t;
//	using configurator_t = typename SERVICE::configurator_t;
//	using input_t = typename SERVICE::input_t;
//	using common_t = typename SERVICE::common_t;
//	using result_t = typename SERVICE::result_t;

	using input_t = typename GenericTypes::input_t;
	using common_t = typename GenericTypes::common_t;
	using result_t = typename GenericTypes::result_t;

	using configurator_t = typename ConfiguratorTypeWrapper_Modes<GenericTypes>::configurator_t_Modes;
	using real_t = typename configurator_t::real_t;

	using service_t = Service<GenericTypes>;
public:
	scATTRACTModes();
	virtual ~scATTRACTModes() {}

	void init(CmdArgs const& args) override;
	void finalize() override;
	void run() override;

private:

	std::shared_ptr<configurator_t> _config;

};

}  // namespace as



#endif /* SRC_SCATTRACT_H_ */
