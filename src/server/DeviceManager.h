/*
 * DeviceManager.h
 *
 *  Created on: May 29, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEMANAGER_H_
#define SRC_DEVICEMANAGER_H_

#ifdef CUDA
#include <memory>
#include <map>
#include <vector>

#include "publicTypes.h"
#include "DeviceManagerHelpers.h"
#include "DeviceItem.h"

namespace as {

class DataItem;
//class DeviceItem;

class DeviceManager {
public:

	DeviceManager() {}
	~DeviceManager() {
		detachAll();
	};

	void attachToDevice(std::shared_ptr<DataItem> item, id_t id, deviceId_t deviceId);
	void detachFromDevice(id_t id, deviceId_t deviceId);
	void detachFromAllDevices(id_t id);
	void detachAllFromDevice(deviceId_t deviceId);
	void detachAll();

	/** returns all deviceIds the item with id is attached to */
	std::vector<deviceId_t> getDeviceIds(id_t id) noexcept;

	/** returns all ids of items that are attached to the device with deviceId */
	std::vector<id_t> getIds(deviceId_t deviceId) noexcept;

	std::shared_ptr<DeviceItem> getItem(id_t id, deviceId_t deviceId) const;

private:
	/** maps deviceId to a set of DataItems (id_t) */
	std::map<deviceId_t, DeviceOccupancy> _deviceOccupancyMap;

	/** maps an id to a set of devices (deviceId_t) the DataItem is attached to */
	DeviceLocationMap _deviceLocationMap;

	std::map<DeviceItemKey, std::shared_ptr<DeviceItem>> _deviceItems;

	/**	checks if item with id is attach to deviceId */
	bool attached(id_t id, deviceId_t deviceId) const;

};

}

#endif // CUDA

#endif /* SRC_DEVICEMANAGER_H_ */
