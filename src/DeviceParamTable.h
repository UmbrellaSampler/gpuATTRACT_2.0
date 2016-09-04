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
		using type_t = typename ParamTable<REAL>::type_t;
		unsigned numTypes;  								/** number of particle/atom types */
		PotShape shape;			/** potential shape that is supported by the table */
		type_t* paramTable;

		type_t getParams(const int& typeA, const int& typeB) const noexcept {
			return paramTable[numTypes*typeA + typeB];
		}
	};

	using HostResc = Desc;

	Desc desc;
	HostResc hostResc;
};

}  // namespace as

#endif

#endif /* SRC_DEVICEPARAMTABLE_H_ */
