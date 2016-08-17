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

template<typename REAL>
class Configurator_6D;

template<typename REAL>
class scATTRACT : public App<REAL> {
public:
	scATTRACT();
	virtual ~scATTRACT() {}

	void init(CmdArgs const& args) override;
	void finalize() override;
	void run() override;

private:
	std::shared_ptr<Configurator_6D<REAL>> _config;

};

}  // namespace as



#endif /* SRC_SCATTRACT_H_ */
