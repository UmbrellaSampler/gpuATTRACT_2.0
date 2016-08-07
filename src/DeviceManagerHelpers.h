/*
 * DeviceManager.tpp
 *
 *  Created on: May 30, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEMANAGERHELPERS_H_
#define SRC_DEVICEMANAGERHELPERS_H_

#ifdef CUDA
#include <stdexcept>
#include <set>
#include <map>
#include "publicTypes.h"

namespace as {

/**
 * Keeps track of which DataItems are attached to a specific device
 */
class DeviceOccupancy {
public:
	void addId(id_t id) noexcept {
		_ids.insert(id);
	}

	bool contains(id_t id) const noexcept {
		return _ids.find(id) != _ids.end();
	}

	std::set<id_t> getIds() const noexcept {
		return _ids;
	}

	void removeId(id_t id) noexcept {
		_ids.erase(id);
	}

private:
	/** IDs of DataItems that are located on a device */
	std::set<id_t> _ids;

};

/**
 * Maps dataItems to a set of deviceIds they are attached to
 */
class DeviceLocationMap {
public:
	void addIdToDevice(id_t id, deviceId_t dId) {
		auto& set = _deviceMap[id]; // create set if not there

//		if (set.find(dId) != set.end()) {
//			std::stringstream what;
//			what << "Cannot attach DataItem. DataItem with id " << id << " is already attached to device " << dId << "\n";
//			throw std::invalid_argument(what.str());
//		}
		set.insert(dId);
	}

	void removeIdFromDevice(id_t id, deviceId_t dId) {
		auto& set = _deviceMap.at(id);
//		if (set.find(dId) == set.end()) {
//			std::stringstream what;
//			what << "Cannot detach DataItem. DataItem with id " << id << " is not attached to device " << dId << "\n";
//			throw std::invalid_argument(what.str());
//		}
		set.erase(set.find(dId));
	}

	std::set<deviceId_t> deviceIds(id_t id) const {
		return _deviceMap.at(id);
	}


private:
	std::map<id_t, std::set<deviceId_t>> _deviceMap;

};

class DeviceItemKey {
public:

	DeviceItemKey(id_t id, deviceId_t deviceId) : _id(id), _deviceId(deviceId) {}

	bool operator == (DeviceItemKey const& rhs) const noexcept {
		return _id == rhs._id && _deviceId == rhs._deviceId;
	}

	bool operator < (DeviceItemKey const& rhs) const noexcept {
		auto p1 = std::pair<id_t, deviceId_t>(_id, _deviceId);
		auto p2 = std::pair<id_t, deviceId_t>(rhs._id, rhs._deviceId);
		return p1 < p2;
	}

private:
	id_t _id;
	deviceId_t _deviceId;
};


} // namespace as
#endif //CUDA



#endif /* SRC_DEVICEMANAGERHELPERS_H_ */
