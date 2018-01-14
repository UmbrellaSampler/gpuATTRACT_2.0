/*
 * readDof_test.cpp
 *
 *  Created on: Jan 13, 2018
 *      Author: uwe
 */


#include <gtest/gtest.h>
#include "readFile.h"
#include "IOException.h"
#include <vector>

using namespace std;
using namespace as;

TEST(readFile, readDof_content) {
	const std::string filename ="./test/resources/readDof_test.dat";
	vector<vector<DOF>> dofs = readDOF(filename);
	ASSERT_EQ(3, dofs.size());

	for(auto const& dof_vec : dofs) {
		ASSERT_EQ(3, dof_vec.size());
		for (auto const& dof : dof_vec) {
			ASSERT_DOUBLE_EQ(1.0, dof._6D.ang.x);
			ASSERT_DOUBLE_EQ(2.0, dof._6D.ang.y);
			ASSERT_DOUBLE_EQ(3.0, dof._6D.ang.z);
			ASSERT_DOUBLE_EQ(4.0, dof._6D.pos.x);
			ASSERT_DOUBLE_EQ(5.0, dof._6D.pos.y);
			ASSERT_DOUBLE_EQ(6.0, dof._6D.pos.z);

			ASSERT_EQ(8, dof.numDofs);
			ASSERT_DOUBLE_EQ(7.0, dof.dofs[0]);
			ASSERT_DOUBLE_EQ(8.0, dof.dofs[1]);
		}
	}
}

TEST(readFile, readDof_wrong_filename) {
	const std::string filename ="non-existing filename";

	try {
		vector<vector<DOF>> dofs = readDOF(filename);
		FAIL();
	} catch (IOException& e) {
		ASSERT_STREQ( std::string("Error reading file '" + filename + "'. The file does not exist").c_str(), e.what() );
	} catch (std::exception& e) {
		FAIL();
	}
}



