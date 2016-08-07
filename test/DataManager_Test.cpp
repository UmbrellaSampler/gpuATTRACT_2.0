/*
 * DataManager_Test.cpp
 *
 *  Created on: May 26, 2016
 *      Author: uwe
 */

#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "Protein.tpp"
#include "DataItem.h"
#include "DataManager.h"


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

