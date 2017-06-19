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

namespace as {

template<typename SERVICE>
class scATTRACT : public App {
	using real_t = typename SERVICE::real_t;
	using configurator_t = typename SERVICE::configurator_t;
	using dof_t = typename SERVICE::dof_t;
	using common_t = typename SERVICE::common_t;
	using result_t = typename SERVICE::result_t;
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
