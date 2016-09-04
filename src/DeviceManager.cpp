/*
 * DeviceManager.cpp
 *
 *  Created on: Jul 19, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include <sstream>
#include <stdexcept>

#include "DeviceManager.h"
#include "DataItem.h"
#include "DeviceItem.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "SimParam.h"
#include "DeviceProtein.h"
#include "DeviceGridUnion.h"
#include "DeviceParamTable.h"
#include "DeviceDataConfigurator.h"
#include "macros.h"
#include "cuda_runtime.h"

using namespace as;

void DeviceManager::attachToDevice(std::shared_ptr<DataItem> item, id_t id, deviceId_t deviceId) {

	if(attached(id, deviceId)) {
		std::stringstream what;
		what << "Attaching failed. DataItem with id " << id << " is already attached at device " << deviceId;
		throw std::invalid_argument(what.str());
	}

	std::shared_ptr<DeviceItem> deviceItem;
	if(auto derived = std::dynamic_pointer_cast<Protein<float>>(item)) {
		deviceItem = std::static_pointer_cast<DeviceItem>(DeviceDataConfigurator::attach(derived, deviceId));
	} else if (auto derived = std::dynamic_pointer_cast<Protein<double>>(item)) {
		deviceItem = std::static_pointer_cast<DeviceItem>(DeviceDataConfigurator::attach(derived, deviceId));
	} else if (auto derived =  std::dynamic_pointer_cast<GridUnion<float>>(item)) {
		deviceItem = std::static_pointer_cast<DeviceItem>(DeviceDataConfigurator::attach(derived, deviceId));
	} else if (auto derived =  std::dynamic_pointer_cast<GridUnion<double>>(item)) {
		deviceItem = std::static_pointer_cast<DeviceItem>(DeviceDataConfigurator::attach(derived, deviceId));
	} else if (auto derived =  std::dynamic_pointer_cast<ParamTable<float>>(item)) {
		deviceItem = std::static_pointer_cast<DeviceItem>(DeviceDataConfigurator::attach(derived, deviceId));
	} else if (auto derived =  std::dynamic_pointer_cast<ParamTable<double>>(item)) {
		deviceItem = std::static_pointer_cast<DeviceItem>(DeviceDataConfigurator::attach(derived, deviceId));
	} else if (std::dynamic_pointer_cast<SimParam<float>>(item)
			|| std::dynamic_pointer_cast<SimParam<double>>(item)) {
		return;
	} else {
		throw std::runtime_error("Invalid DataItem. Dynamic cast failed");
	}
	_deviceOccupancyMap[deviceId].addId(id);
	_deviceLocationMap.addIdToDevice(id, deviceId);
	_deviceItems[DeviceItemKey(id, deviceId)] = deviceItem;
}

void DeviceManager::detachFromDevice(id_t id, deviceId_t deviceId) {

	if(!attached(id, deviceId)) {
		std::stringstream what;
		what << "Detaching failed. DataItem with id " << id << " is not attached at device " << deviceId;
		throw std::invalid_argument(what.str());
	}

	auto& item = _deviceItems.at(DeviceItemKey(id, deviceId));

	if(auto derived = std::dynamic_pointer_cast<DeviceProtein<float>>(item)) {
		DeviceDataConfigurator::detach(derived, deviceId);
	} else if (auto derived = std::dynamic_pointer_cast<DeviceProtein<double>>(item)) {
		DeviceDataConfigurator::detach(derived, deviceId);
	} else if (auto derived = std::dynamic_pointer_cast<DeviceGridUnion<float>>(item)) {
		DeviceDataConfigurator::detach(derived, deviceId);
	} else if (auto derived = std::dynamic_pointer_cast<DeviceGridUnion<double>>(item)) {
		DeviceDataConfigurator::detach(derived, deviceId);
	} else if (auto derived = std::dynamic_pointer_cast<DeviceParamTable<float>>(item)) {
		DeviceDataConfigurator::detach(derived, deviceId);
	} else if (auto derived = std::dynamic_pointer_cast<DeviceParamTable<double>>(item)) {
		DeviceDataConfigurator::detach(derived, deviceId);
	} else {
		throw std::runtime_error("Invalid DeviceItem. Dynamic cast failed");
	}

	_deviceOccupancyMap[deviceId].removeId(id);
	_deviceLocationMap.removeIdFromDevice(id, deviceId);
	_deviceItems.erase(DeviceItemKey(id, deviceId));
}

void DeviceManager::detachFromAllDevices(id_t id) {
	const auto deviceIds = _deviceLocationMap.deviceIds(id);
	for (auto deviceId : deviceIds) {
		detachFromDevice(id, deviceId);
	}
}

void DeviceManager::detachAllFromDevice(deviceId_t deviceId) {
	const auto ids = _deviceOccupancyMap[deviceId].getIds();
	for (auto id : ids) {
		detachFromDevice(id, deviceId);
	}
}

void DeviceManager::detachAll() {
	const auto numDevices = _deviceOccupancyMap.size();
	for (size_t deviceId = 0; deviceId < numDevices; ++deviceId) {
		detachAllFromDevice(deviceId);
	}
}

bool DeviceManager::attached(id_t id, deviceId_t deviceId) const {
	return _deviceItems.find(DeviceItemKey(id, deviceId)) != _deviceItems.end();
}

std::vector<deviceId_t> DeviceManager::getDeviceIds(id_t id) noexcept {
	const auto deviceIdsSet = _deviceLocationMap.deviceIds(id);
	std::vector<deviceId_t> deviceIds(deviceIdsSet.begin(), deviceIdsSet.end());
	return deviceIds;
}

std::vector<id_t> DeviceManager::getIds(deviceId_t deviceId) noexcept {
	const auto idsSet = _deviceOccupancyMap[deviceId].getIds();
	std::vector<id_t> ids(idsSet.begin(), idsSet.end());
	return ids;
}

std::shared_ptr<DeviceItem> DeviceManager::getItem(id_t id, deviceId_t deviceId) const {
	if(!attached(id, deviceId)) {
		std::stringstream what;
		what << "Call to getItem() failed. DataItem with id " << id << " is not attached at device " << deviceId;
		throw std::invalid_argument(what.str());
	}
	return _deviceItems.at(DeviceItemKey(id, deviceId)); // does not throw due to the check before
}


#endif
