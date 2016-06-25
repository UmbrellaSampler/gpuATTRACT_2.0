/*
 * DeviceManager.h
 *
 *  Created on: May 29, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEMANAGER_H_
#define SRC_DEVICEMANAGER_H_

#include "DeviceProtein.h"

namespace as {

template <typename REAL>
class Protein;

template <typename REAL>
class DeviceManager {
public:
	static DeviceProtein<REAL> createDeviceProtein(Protein<REAL> const* protein);
};

}

#endif /* SRC_DEVICEMANAGER_H_ */
