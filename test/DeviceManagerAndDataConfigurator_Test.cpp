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
#include "readFile.h"
#include "DeviceDataConfigurator.h"
#include "TypeMap.h"

using namespace as;

class DeviceManagerAndDataConfig_TEST : public ::testing::Test {
protected: // could also be public according to gtest Primer.md
	virtual void SetUp() {
		_items = {
				getTestProtein<float>(),
				getTestProtein<double>(),
				getTestGrid<float>(),
				getTestGrid<double>()
			};
		_ids = {0, 1, 2, 3};

		auto prot = _items[0];

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

#endif



