/*
 * CPUEnergyService6DTest.cpp
 *
 *  Created on: Jan 21, 2018
 *      Author: uwe
 */

#include <gtest/gtest.h>
#include <string>
#include <iostream>

#include "Protein.h"
#include "GridUnion.h"
#include "SimParam.h"
#include "ParamTable.h"
#include "CmdParser.h"
#include "CmdArgs.h"
#include "StringToArgConverter.h"
#include "Types_6D.h"
#include "Service.h"
#include "ConfiguratorTypeWrapper.h"
#include "Request.h"
#include "Server.h"
#include "readFile.h"


using namespace as;
using namespace std;
//using ::testing::;

TEST(Scoring, CPU) {
	constexpr double eps = 0.0005;

	using GenericTypes = Types_6D<double>;
	using input_t = typename GenericTypes::input_t;
	using common_t = typename GenericTypes::common_t;
	using result_t = typename GenericTypes::result_t;

	using configurator_t = typename ConfiguratorTypeWrapper<GenericTypes>::configurator_t;
	const string path = "./test/resources/";

	const string receptorPdb = path + "receptorr.pdb";
	const string ligandPdb = path + "ligandr.pdb";
	const string gridFile = path + "receptorgrid_test.grid";
	const string gridAlphabet = path + "receptorgrid.alphabet";
	const string paramFile  = path + "attract.par";
	const string dofFileSystsearch = path + "systsearch.dat";
	const string dofFileDocking = path + "out_docking_orig.dat";
	const string scoreFileSystsearch = path + "out_systsearch_orig.score";
	const string scoreFileDocking = path + "out_docking_orig.score";
//	const string scoreFileSystsearch = path + "out_attract_modes_no_fixRec.score";

	const string cmd = string("command")
			+ " sc "
			+ " --dof " 		+ dofFileDocking
			+ " -r " 			+ receptorPdb
			+ " -l " 			+ ligandPdb
			+ " --gridrec " 	+ gridFile
			+ " -p " 			+ paramFile
			+ " --alphabetrec " + gridAlphabet
			+ " -c " 			+ "1"
			+ " --prec " 		+ "double";

	StringToArgConverter converter(cmd);

	CmdParser parser;
	parser.parse(converter.argc(), converter.argv());

	CmdArgs cmdArgs = parser.args();

	configurator_t serverConfigurator;
	serverConfigurator.init(cmdArgs);

	auto& dofs = serverConfigurator.dofs();
	auto server = serverConfigurator.server();
	auto& common = serverConfigurator.common();
	size_t numDofs = dofs.size();
	Request<input_t, common_t> request(dofs.data(), numDofs, common);
	std::shared_ptr<Request<input_t, common_t>> shared_request= std::make_shared<Request<input_t,common_t>>(request);
	server->submit(shared_request);

	auto results = std::vector<result_t>(dofs.size());
	server->wait(shared_request, results.data());

	vector<Result> refResults = readResult(scoreFileDocking);

	ASSERT_EQ(refResults.size(), results.size());
	bool passed = true;
	for (size_t i = 0; i < refResults.size(); ++i) {
		if (abs((refResults[i].E -results[i].E)/refResults[i].E) > eps) {
			passed = false;
			cout << "s" << i << " E:  " << refResults[i].E << " " << results[i].E << endl;
			cout << "s" << i << " G0: " << refResults[i].gradients[0]._6D[0] << " " << results[i].ang.x << endl;
			cout << "s" << i << " G1: " << refResults[i].gradients[0]._6D[1] << " " << results[i].ang.y << endl;
			cout << "s" << i << " G2: " << refResults[i].gradients[0]._6D[2] << " " << results[i].ang.z << endl;
			cout << "s" << i << " G3: " << refResults[i].gradients[0]._6D[3] << " " << results[i].pos.x << endl;
			cout << "s" << i << " G4: " << refResults[i].gradients[0]._6D[4] << " " << results[i].pos.y << endl;
			cout << "s" << i << " G5: " << refResults[i].gradients[0]._6D[5] << " " << results[i].pos.z << endl;
			cout << endl;
		}

	}
	ASSERT_TRUE(passed);
}



