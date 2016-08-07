/*
 * DataManager.tpp
 *
 *  Created on: May 22, 2016
 *      Author: uwe
 */

#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cassert>
#include "DataManager.h"
#include "DataItem.h"

#ifdef CUDA
#include "DeviceManager.h"
#include "cuda_runtime.h"
#endif


namespace as {

id_t DataManager::add(std::shared_ptr<DataItem> dataItem) {
	if (dataItem.get() == nullptr) {
		throw std::invalid_argument("Invalid DataItem (nullptr).");
	}
	_dataItems.push_back(dataItem);
	return _dataItems.size()-1;
}

std::vector<id_t> DataManager::add(std::vector<std::shared_ptr<DataItem>> const& cont) {
	std::vector<id_t> ids;
	for (auto& item : cont) {
		auto id = add(item);
		ids.push_back(id);
	}
	return ids;
}

bool DataManager::isValid(id_t id) {
	return (id < _dataItems.size()) && (_dataItems[id] != nullptr);
}

std::shared_ptr<DataItem> DataManager::get(id_t id) {
	if (!isValid(id)) {
		std::stringstream stream;
		stream << id;
		throw std::invalid_argument("Invalid id for DataItem (id = "+ stream.str() + ").");
	}
	return _dataItems[id];
}

void DataManager::remove(id_t id) {
	if (!isValid(id)) {
		std::stringstream stream;
		stream << id;
		throw std::invalid_argument("Invalid id for DataItem (id = "+ stream.str() + ").");
	}
	_dataItems[id] = nullptr;
}

void DataManager::remove(std::vector<id_t> const& cont) {
	for (auto& id : cont) {
		remove(id);
	}
}

void DataManager::removeAll() {
	_dataItems.clear();
}


#ifdef CUDA
void DataManager::attachToDevice(id_t id, deviceId_t deviceId) {
	auto& item = _dataItems[id];
	_deviceManager.attachToDevice(item, id, deviceId);
}
void DataManager::detachFromDevice(id_t id, deviceId_t deviceId) {
	_deviceManager.detachFromDevice(id, deviceId);
}

void DataManager::attachAllDataToDevice(deviceId_t deviceId) {
	for (size_t id = 0; id < _dataItems.size(); ++id) {
		if (isValid(id)) {
			attachToDevice(id, deviceId);
		}
	}
}

void DataManager::releaseDevice(deviceId_t deviceId) {
	_deviceManager.detachAllFromDevice(deviceId);
}
void DataManager::releaseAllDevices() {
	_deviceManager.detachAll();
}

std::vector<deviceId_t> DataManager::getCommonDeviceIds(std::vector<id_t> const& ids) {
	std::vector<std::vector<deviceId_t>> idVecs;
	for (size_t i = 0; i < ids.size(); ++i) {
		auto idsOfItem = _deviceManager.getDeviceIds(ids[i]);
		assert(std::is_sorted(idsOfItem.begin(), idsOfItem.end()));
		idVecs.push_back(std::move(idsOfItem));
	}
	/* build intersection */
	/* initialize with deviceIds of first item and build intersection with remaining */
	std::set<deviceId_t> intersect(idVecs[0].begin(), idVecs[0].end());
	for (size_t i = 1; i < idVecs.size(); ++i) {
		auto const& idVec = idVecs[i];
		std::set<deviceId_t> tmp;
		std::set_intersection(intersect.begin(), intersect.end(), idVec.begin(), idVec.end(),
				std::insert_iterator<std::set<deviceId_t>>(tmp, tmp.begin()));
		intersect = std::move(tmp);
	}

	return std::vector<deviceId_t>(intersect.begin(), intersect.end());
}

#endif //CUDA

} // namespace as

