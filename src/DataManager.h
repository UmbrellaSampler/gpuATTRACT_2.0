/*
 * DataManager.h
 *
 *  Created on: Mar 30, 2016
 *      Author: uwe
 */

#ifndef SRC_DATAMANAGER_H_
#define SRC_DATAMANAGER_H_

#include <set>
#include <map>
#include <vector>
#include <memory>

#include "publicTypes.h"
#include "DeviceManager.h"

namespace as {

class DataItem;

class DataManager {

public:
	DataManager() {};
	~DataManager() {};

	id_t add(std::shared_ptr<DataItem> dataItem);

	std::vector<id_t> add(std::vector<std::shared_ptr<DataItem>> const& cont);

	std::shared_ptr<DataItem> get(id_t);

	void remove(id_t id);
	void remove(std::vector<id_t> const& cont);

	/**
	 * Removes all DataItems. Using old/previous ids on the CPU leads now to undefined behavior;
	 * You can still use them to address items on the gpu, that are already attached.
	 */
	void removeAll();

private:
	bool isValid(id_t id);
	std::vector<std::shared_ptr<DataItem>> _dataItems;

#ifdef CUDA
public:
	void attachToDevice(id_t id, deviceId_t deviceId);
	void detachFromDevice(id_t id, deviceId_t deviceId);

	void attachAllDataToDevice(deviceId_t deviceId);
	void attachDataToAllDevices();

	void releaseDevice(deviceId_t deviceId);
	void releaseAllDevices();

	std::vector<deviceId_t> getCommonDeviceIds(std::vector<id_t> const& ids);

	DeviceManager _deviceManager;

#endif

};

}  // namespace as



#endif /* SRC_DATAMANAGER_H_ */
