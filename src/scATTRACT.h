/*
 * scATTRACT.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#ifndef SRC_SCATTRACT_H_
#define SRC_SCATTRACT_H_

#include "App.h"

namespace as {

template<typename REAL>
class scATTRACT : public App<REAL> {

	void run(CmdArgs const& args);
};

}  // namespace as



#endif /* SRC_SCATTRACT_H_ */
