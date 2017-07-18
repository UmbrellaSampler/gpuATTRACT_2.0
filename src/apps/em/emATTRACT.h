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

namespace as {

template<typename SERVICE>
class emATTRACT : public App {
	using real_t = typename SERVICE::real_t;
	using configurator_t = typename SERVICE::configurator_t;
	using dof_t = typename SERVICE::dof_t;
	using common_t = typename SERVICE::common_t;
	using result_t = typename SERVICE::result_t;
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
