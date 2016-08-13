/*
 * DeviceParamTable.h
 *
 *  Created on: Aug 7, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEPARAMTABLE_H_
#define SRC_DEVICEPARAMTABLE_H_

#ifdef CUDA

#include "DeviceItem.h"
#include "ParamTable.h"

namespace as {

template <typename REAL>
class DeviceParamTable : public DeviceItem {
public:
	struct Desc {
		unsigned numTypes;  								/** number of particle/atom types */
		PotShape shape;			/** potential shape that is supported by the table */
		typename ParamTable<REAL>::type_t* paramTable;
	};

	using HostResc = Desc;

	Desc desc;
	HostResc hostResc;
};

}  // namespace as

#endif

#endif /* SRC_DEVICEPARAMTABLE_H_ */
