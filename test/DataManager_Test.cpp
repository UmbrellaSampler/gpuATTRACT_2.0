/*
 * DataManager_Test.cpp
 *
 *  Created on: May 26, 2016
 *      Author: uwe
 */

#include <gtest/gtest.h>
#include <memory>
#include "DataManager.tpp"
#include "Protein.tpp"


TEST(DataManager, addProtein) {
	using namespace as;
	DataManager<float> dataMng;
	std::shared_ptr<Protein<float>> protein;
	dataMng.addProtein(protein);
	try {
		dataMng.addProtein(protein);
		FAIL();
	} catch (std::invalid_argument& e) {

		EXPECT_TRUE(true);
	} catch (std::exception& e) {
		// invalid_argument exception expected
		FAIL();
	}

	for (int i = 0; i < 10; ++i) {
		std::shared_ptr<Protein<float>> protein = std::make_shared<Protein<float>>();
		proteinId_t id = dataMng.addProtein(protein);
		EXPECT_EQ(id, i);
	}
}

