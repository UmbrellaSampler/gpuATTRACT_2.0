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

#ifndef __CUDACC__

#ifndef __device__
#define __device__
#endif

#endif

namespace as {

template <typename REAL>
class DeviceParamTable : public DeviceItem {
public:
	struct Desc {
		using type_t = typename ParamTable<REAL>::type_t;
		unsigned numTypes;  								/** number of particle/atom types */
		PotShape shape;			/** potential shape that is supported by the table */
		type_t* paramTable;

		__device__
		type_t getParams(const int& typeA, const int& typeB) const noexcept {
			return paramTable[numTypes*typeA + typeB];
		}
	};

	using HostResc = Desc;

	Desc desc;
	HostResc hostResc;
};

template<typename REAL>
using d_ParamTable = typename DeviceParamTable<REAL>::Desc;

}  // namespace as

#endif

#endif /* SRC_DEVICEPARAMTABLE_H_ */
