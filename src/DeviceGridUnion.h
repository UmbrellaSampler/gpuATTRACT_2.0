/*
 * DeviceGridUnion.h
 *
 *  Created on: May 29, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEGRIDUNION_H_
#define SRC_DEVICEGRIDUNION_H_

#ifdef CUDA
#include <memory>
#include "DeviceItem.h"
#include "DeviceIntrplGrid.h"
#include "DeviceNLGrid.h"
/**
 * DeviceProtein represents a Protein on a device.
 */
namespace as {

template <typename REAL>
class DeviceGridUnion : public DeviceItem {
public:

	struct Desc {
		typename DeviceIntrplGrid<REAL>::Desc inner;
		typename DeviceIntrplGrid<REAL>::Desc outer;
		typename DeviceNLGrid<REAL>::Desc NL;
	};

	struct HostResc {
		typename DeviceIntrplGrid<REAL>::HostResc inner;
		typename DeviceIntrplGrid<REAL>::HostResc outer;
		typename DeviceNLGrid<REAL>::HostResc NL;
	};


	Desc getDesc() {
		Desc desc;
		desc.inner = inner->desc;
		desc.outer = outer->desc;
		desc.NL = NL->desc;
		return desc;
	}

	HostResc getHostResc() {
		HostResc resc;
		resc.inner = inner->hostResc;
		resc.outer = outer->hostResc;
		resc.NL = NL->hostResc;
		return resc;
	}

	std::shared_ptr<DeviceIntrplGrid<REAL>> inner;
	std::shared_ptr<DeviceIntrplGrid<REAL>> outer;
	std::shared_ptr<DeviceNLGrid<REAL>> NL;

};


}

#endif // CUDA



#endif /* SRC_DEVICEGRIDUNION_H_ */
