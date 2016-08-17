/*
 * scATTRACT.tpp
 *
 *  Created on: Aug 17, 2016
 *      Author: uwe
 */

#ifndef SRC_SCATTRACT_TPP_
#define SRC_SCATTRACT_TPP_

#include "scATTRACT.h"
#include "Configurator_6D.h"

namespace as {

template<typename REAL>
scATTRACT<REAL>::scATTRACT() : _config(new Configurator_6D<REAL>()) {}

template<typename REAL>
void scATTRACT<REAL>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename REAL>
void scATTRACT<REAL>::finalize() {
	_config->finalize();
}

template<typename REAL>
void scATTRACT<REAL>::run() {

}

}  // namespace as



#endif /* SRC_SCATTRACT_TPP_ */
