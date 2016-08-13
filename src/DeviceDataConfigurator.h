/*
 * CUDADataConfigurator.cpp
 *
 *  Created on: Jul 19, 2016
 *      Author: uwe
 */

#ifndef SRC_CUDADATACONFIGURATOR_H_
#define SRC_CUDADATACONFIGURATOR_H_

#ifdef CUDA

#include <memory>
#include "publicTypes.h"


namespace as {

template<typename REAL>
class DeviceProtein;

template<typename REAL>
class Protein;

template<typename REAL>
class DeviceGridUnion;

template<typename REAL>
class GridUnion;

template<typename REAL>
class IntrplGrid;

template<typename REAL>
class DeviceIntrplGrid;

template<typename REAL>
class NLGrid;

template<typename REAL>
class DeviceNLGrid;

template<typename REAL>
class ParamTable;

template<typename REAL>
class DeviceParamTable;

class DeviceDataConfigurator {

public:
	template<typename REAL>
	static std::shared_ptr<DeviceProtein<REAL>> attach(const std::shared_ptr<Protein<REAL>>, deviceId_t);

	template<typename REAL>
	static void detach(const std::shared_ptr<DeviceProtein<REAL>>, deviceId_t);

	template<typename REAL>
	static std::shared_ptr<DeviceGridUnion<REAL>> attach(const std::shared_ptr<GridUnion<REAL>> , deviceId_t);

	template<typename REAL>
	static void detach(const std::shared_ptr<DeviceGridUnion<REAL>> , deviceId_t);

	template<typename REAL>
	static std::shared_ptr<DeviceParamTable<REAL>> attach(const std::shared_ptr<ParamTable<REAL>>, deviceId_t);

	template<typename REAL>
	static void detach(const std::shared_ptr<DeviceParamTable<REAL>>, deviceId_t);

	static void setEnableDeviceCheck(bool status) {
		enableDeviceCheck = status;
	}

private:
	template<typename REAL>
	static std::shared_ptr<DeviceIntrplGrid<REAL>> attach(const std::shared_ptr<IntrplGrid<REAL>> protein);

	template<typename REAL>
	static void detach(const std::shared_ptr<DeviceIntrplGrid<REAL>> protein);

	template<typename REAL>
	static std::shared_ptr<DeviceNLGrid<REAL>> attach(const std::shared_ptr<NLGrid<REAL>> protein);

	template<typename REAL>
	static void detach(const std::shared_ptr<DeviceNLGrid<REAL>> protein);

	/**
	 * Checks if a cuda device with deviceId is available. If yes, this device is
	 * set to the current device. Otherwise, an exception is thrown.
	 */
	static void checkDeviceIdAndSetCurrent(deviceId_t deviceId);


	/** For testing purposes. Must be set to true under release conditions, which is the default **/
	static bool enableDeviceCheck;
};


}  // namespace as

#endif

#endif /* SRC_CUDADATACONFIGURATOR_CPP_ */
