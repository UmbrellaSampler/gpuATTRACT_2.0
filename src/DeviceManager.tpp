/*
 * DeviceManager.tpp
 *
 *  Created on: May 30, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEMANAGER_TPP_
#define SRC_DEVICEMANAGER_TPP_



#include "DeviceManager.h"

#ifdef CUDA
using namespace as;

template <typename REAL>
DeviceProtein<REAL> DeviceManager<REAL>::createDeviceProtein(Protein<REAL> const* protein) {
	DeviceProtein<REAL> deviceProtein;
	return deviceProtein;

}

#endif


#endif /* SRC_DEVICEMANAGER_TPP_ */
