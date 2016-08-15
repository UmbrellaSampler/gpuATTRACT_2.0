/*
 * DeviceManager_Test.cpp
 *
 *  Created on: May 30, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include <gtest/gtest.h>
#include <memory>
#include "DeviceManager.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "../src/fileIO/readFile.h"
#include "DeviceDataConfigurator.h"

using namespace as;

class DeviceManagerAndDataConfig_TEST : public ::testing::Test {
protected: // could also be public according to gtest Primer.md
	virtual void SetUp() {
		_items = {
				getTestProtein<float>(),
				getTestProtein<double>(),
				getTestGrid<float>(),
				getTestGrid<double>(),
				getTestTable<float>(),
				getTestTable<double>()
			};
		_ids = std::vector<id_t>(_items.size());
		std::iota(_ids.begin(), _ids.end(), 0);
	}

	virtual void TearDown() {}

	template<typename REAL>
	std::shared_ptr<Protein<REAL>> getTestProtein() {
		auto protein = createProteinFromPDB<REAL>(pdbName);
		protein->setNumMappedTypes(1);
		protein->getOrCreateMappedPtr();
		applyDefaultMapping(protein->numAtoms(), protein->type(), protein->mappedType());
		return protein;
	}

	template<typename REAL>
	std::shared_ptr<GridUnion<REAL>> getTestGrid() {
		return createGridFromGridFile<REAL>(gridName);
	}

	template<typename REAL>
	std::shared_ptr<ParamTable<REAL>> getTestTable() {
		return createParamTableFromFile<REAL>(tableName);
	}

	std::vector<std::shared_ptr<DataItem>>& getItems_ref() {
		return _items;
	}

	std::vector<id_t>& getIds_ref() {
		return _ids;
	}
private:

	std::vector<std::shared_ptr<DataItem>> _items;
	std::vector<id_t> _ids;
	const std::string pdbName = "/home/uwe/eclipse/gpuATTRACT/gpuATTRACT_2.0/test/resources/ligandr.pdb";
	const std::string gridName = "/home/uwe/eclipse/gpuATTRACT/gpuATTRACT_2.0/test/resources/receptorgrid.grid";
	const std::string tableName = "/home/uwe/eclipse/gpuATTRACT/gpuATTRACT_2.0/test/resources/attract.par";


};

TEST_F(DeviceManagerAndDataConfig_TEST, DeviceDataConfigurator_attach_detach) {

	const deviceId_t deviceId = 0;
	auto& items = getItems_ref();

	try {
		auto deviceData0 = DeviceDataConfigurator::attach(std::dynamic_pointer_cast<Protein<float>>(items[0]), deviceId);
		DeviceDataConfigurator::detach(deviceData0, deviceId);
		auto deviceData1 = DeviceDataConfigurator::attach(std::dynamic_pointer_cast<Protein<double>>(items[1]), deviceId);
		DeviceDataConfigurator::detach(deviceData1, deviceId);
		auto deviceData2 = DeviceDataConfigurator::attach(std::dynamic_pointer_cast<GridUnion<float>>(items[2]), deviceId);
		DeviceDataConfigurator::detach(deviceData2, deviceId);
		auto deviceData3 = DeviceDataConfigurator::attach(std::dynamic_pointer_cast<GridUnion<double>>(items[3]), deviceId);
		DeviceDataConfigurator::detach(deviceData3, deviceId);
	} catch (std::exception& e) {
		FAIL() << e.what();
	}
	SUCCEED();
}


TEST_F(DeviceManagerAndDataConfig_TEST, DeviceManager_attach_detach) {

	DeviceManager mng;
	auto& items = getItems_ref();
	auto& ids = getIds_ref();

	for(size_t i = 0; i < items.size(); ++i) {
		mng.attachToDevice(items[i],ids[i],0);
	}

	for(size_t i = 0; i < items.size(); ++i) {
		try {
			mng.attachToDevice(items[i],ids[i],0);
		} catch (std::invalid_argument& e){
			std::stringstream what;
			what << "Attaching failed. DataItem with id " << ids[i] << " is already attached at device " << 0;
			ASSERT_STREQ( what.str().c_str(), e.what());
		} catch (std::exception& e) {
			FAIL() << e.what();
		}
	}

	for(size_t i = 0; i < items.size(); ++i) {
		mng.detachFromDevice(ids[i],0);
	}

	for(size_t i = 0; i < items.size(); ++i) {
		try {
			mng.detachFromDevice(ids[i],0);
		} catch (std::invalid_argument& e){
			std::stringstream what;
			what << "Detaching failed. DataItem with id " << ids[i] << " is not attached at device " << 0;
			ASSERT_STREQ( what.str().c_str(), e.what());
		} catch (std::exception& e) {
			FAIL() << e.what();
		}
	}

	class Dummy : public DataItem {

	};

	auto dummy = std::make_shared<Dummy>();
	try {
		mng.attachToDevice(dummy,0,0);
	} catch (std::runtime_error& e) {
		ASSERT_STREQ( "Invalid DataItem. Dynamic cast failed", e.what());
	} catch (std::exception& e) {
		FAIL() << e.what();
	}

}

TEST_F(DeviceManagerAndDataConfig_TEST, DeviceManager_getIds) {
	DeviceManager mng;
	auto& items = getItems_ref();
	auto& ids = getIds_ref();

	DeviceDataConfigurator::setEnableDeviceCheck(false);

	const std::vector<std::vector<deviceId_t>> deviceIds = {
			{0,1,2}, // item 0 goes to device 0,1,2
			{2,3,4}, // item 1 goes to device 2,3,4 ...
			{3,4},
			{0,1}
	};

	const std::vector<std::vector<id_t>> commonIds = {
			{ids[0],ids[3]}, // device 0 gets item 0,3
			{ids[0],ids[3]}, // device 1 gets itme 0,3 ...
			{ids[0],ids[1]},
			{ids[1],ids[2]},
			{ids[1],ids[2]}
	};

	for (id_t id = 0; id < deviceIds.size(); ++id) {
		for (auto deviceId : deviceIds[id]) {
			mng.attachToDevice(items[id], ids[id], deviceId);
		}
	}

	for (id_t id = 0; id < deviceIds.size(); ++id) {
		auto deviceIds_ = mng.getDeviceIds(id);
		EXPECT_EQ(deviceIds[id],deviceIds_);
	}

	for (deviceId_t deviceId = 0; deviceId < commonIds.size(); ++deviceId) {
		auto ids = mng.getIds(deviceId); // get all item ids on device with deviceId
		EXPECT_EQ(commonIds[deviceId],ids);
	}

	mng.detachAll();
	DeviceDataConfigurator::setEnableDeviceCheck(true);
}

#endif



