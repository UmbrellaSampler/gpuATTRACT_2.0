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


using namespace as;
using namespace std;
//using ::testing::;

TEST(Scoring, CPU) {

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

//	cout << cmd << endl;

	StringToArgConverter converter(cmd);

//	cout << converter.argc() << endl;

//	for (int i = 0; i < converter.argc(); ++i) {
//		cout << converter.argv()[i] << " ";
//	}
//	cout << endl;

	CmdParser parser;
	parser.parse(converter.argc(), converter.argv());

	CmdArgs cmdArgs = parser.args();

//	cout << appParams << endl;

	configurator_t serverConfigurator;
	serverConfigurator.init(cmdArgs);

	auto& dofs = serverConfigurator.dofs();
	auto server = serverConfigurator.server();
	auto& common = serverConfigurator.common();
	size_t numDofs = dofs.size();
	Request<input_t, common_t> request(dofs.data(), numDofs, common);
	server->submit(request);

	auto results = std::vector<result_t>(dofs.size());
	server->wait(request, results.data());

	std::cout << results[0] << std::endl;





//	shared_ptr<Protein<double>> receptor = createProteinFromPDB<Protein<double>>(receptorPdb);
//	shared_ptr<Protein<double>> ligand = createProteinFromPDB<Protein<double>>(ligandPdb);
//	shared_ptr<GridUnion<double>> grid = createGridFromGridFile<GridUnion<double>>(gridFile);
//	shared
//
//
//	CPUEnergyService6D<double> service(make_shared<DataManager>());
//
//	CPUEnergyService6D<double>::itemProcessor_t itemProcessor = service.createItemProcessor();
//	unique_ptr<CPUEnergyService6D<double>::workItem_t> workItem = make_unique<CPUEnergyService6D<double>::workItem_t>();
//
//	itemProcessor(workItem.get());
}



