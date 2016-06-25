/*
 * DeviceManager_Test.cpp
 *
 *  Created on: May 30, 2016
 *      Author: uwe
 */

#include <gtest/gtest.h>
#include "DeviceManager.tpp"
#include "Protein.tpp"
#include "DeviceProtein.h"

using namespace as;

TEST(DeviceManager, createDeviceProtein) {
	Protein<float> protein;

	auto devProtein = DeviceManager<float>::createDeviceProtein(&protein);
	EXPECT_TRUE(true);
}



