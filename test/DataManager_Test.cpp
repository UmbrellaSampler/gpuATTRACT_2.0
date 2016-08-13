/*
 * DataManager_Test.cpp
 *
 *  Created on: May 26, 2016
 *      Author: uwe
 */

#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "DataItem.h"
#include "DataManager.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "readFile.h"
#include "DeviceDataConfigurator.h"

using namespace as;

class DataManager_TEST : public ::testing::Test {
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

private:

	std::vector<std::shared_ptr<DataItem>> _items;
	const std::string pdbName = "/home/uwe/eclipse/gpuATTRACT/gpuATTRACT_2.0/test/resources/ligandr.pdb";
	const std::string gridName = "/home/uwe/eclipse/gpuATTRACT/gpuATTRACT_2.0/test/resources/receptorgrid.grid";
	const std::string tableName = "/home/uwe/eclipse/gpuATTRACT/gpuATTRACT_2.0/test/resources/attract.par";


};

TEST(DataManager, add) {
	using namespace as;
	DataManager dataMng;
	std::shared_ptr<Protein<float>> protein;

	/*
	 * TEST: nullptr
	 */
	try {
		dataMng.add(protein);
		FAIL();
	} catch (std::invalid_argument& e) {
		ASSERT_STREQ( "Invalid DataItem (nullptr).", e.what() );
	} catch (std::exception& e) {
		// invalid_argument exception expected
		FAIL();
	}

	/*
	 * TEST: add individual valid items
	 */
	int size = 10;
	std::vector<std::shared_ptr<DataItem>> vec;
	for (int i = 0; i < size; ++i) {
		std::shared_ptr<DataItem> protein = std::make_shared<Protein<float>>();
		vec.push_back(protein);
	}

	for (int i = 0; i < size; ++i) {
		id_t id = dataMng.add(vec[i]);
		EXPECT_EQ(id, i);
	}

	/*
	 * TEST: add valid collection
	 */
	auto ids = dataMng.add(vec);
	for (int i = 0; i < size; ++i) {
		EXPECT_EQ(ids[i], size + i);
	}

}

#ifdef CUDA

TEST_F(DataManager_TEST, getCommonDeviceIds) {
	DataManager mng;
	auto& items = getItems_ref();
	auto ids = mng.add(items);

	const std::vector<std::vector<deviceId_t>> deviceIdsVec = {
			{0,1,2}, // item 0 goes to device 0,1,2
			{0,2,3}, // item 1 goes to device 1,2,3 ...
			{0,3},
			{0,1,2,3}
	};

	const std::vector<std::vector<id_t>> itemsVec = {
			{0,1,2,3},
			{0,1,3},
			{2,3},
			{0},
			{1,2,3}
	};

	// expected result
	const std::vector<std::vector<deviceId_t>> commonDevicesVec = {
			{0},
			{0,2},
			{0,3},
			{0,1,2},
			{0,3}
	};

	DeviceDataConfigurator::setEnableDeviceCheck(false);

	for (id_t id = 0; id < deviceIdsVec.size(); ++id) {
		for (auto deviceId : deviceIdsVec[id]) {
			mng.attachToDevice(ids[id], deviceId);
		}
	}

	int j = 0;
	for (auto const& itemIds : itemsVec) {
		auto commonIds = mng.getCommonDeviceIds(itemIds);
		EXPECT_EQ(commonDevicesVec[j++], commonIds);
	}


	mng.releaseAllDevices();
	DeviceDataConfigurator::setEnableDeviceCheck(true);



}

#endif

