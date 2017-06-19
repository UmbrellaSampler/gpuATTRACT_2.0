/*
 * DeviceItem.h
 *
 *  Created on: Jul 19, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEITEM_H_
#define SRC_DEVICEITEM_H_

#ifdef CUDA
#include "publicTypes.h"

namespace as {

class DeviceItem {

protected:
	DeviceItem() {}
	virtual ~DeviceItem() {}
};

}  // namespace as

#endif



#endif /* SRC_DEVICEITEM_H_ */
