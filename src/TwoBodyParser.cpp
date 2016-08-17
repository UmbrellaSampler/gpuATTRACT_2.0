/*
 * TwoBodyParser.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <string>
#include "TwoBodyParser.h"
#include "CmdParserHelper.h"

using namespace as;
using namespace std;
namespace po = boost::program_options;

vector<po::options_description> TwoBodyParser::options() const noexcept {
	std::vector<po::options_description> options;

	po::options_description config("generic");
	config.add_options()
			("help", "print this help massage")
			("config", po::value<string>(), "configuration file")
			("prec", po::value<string>()->default_value("single"), "arithmetic precision ('single', 'double')");
	options.push_back(config);

	po::options_description input("input files");
	input.add_options()
			("dof"     			  , po::value<string>()->required()								, "structure (DOF) file")
			("receptor-pdb,r"     , po::value<string>()->default_value("receptorr.pdb")			, "pdb-file of receptor")
			("ligand-pdb,l"       , po::value<string>()->default_value("ligand.pdb")   			, "pdb-file of ligand")
			("grid,g"             , po::value<string>()->default_value("receptorgrid.grid")		, "receptor grid file")
			("par,p"	          , po::value<string>()->default_value("attract.par")			, "attract forcefield parameter file")
			("alphabet,a"		  , po::value<string>()->default_value("receptorgrid.alphabet")	, "receptor grid alphabet file");
	options.push_back(input);

	po::options_description concurrency("concurrency");
#ifndef CUDA
	concurrency.add_options()
			("numCPUs,c", po::value<int>()->default_value(1), "number of CPU threads for CPU mode")
			("chunkSize", po::value<int>()->default_value(1000), "number of concurrently processed structures at the server");
#else
	concurrency.add_options()
			("numCPUs,c", po::value<int>()->default_value(0), "number of CPU threads for CPU mode")
			("device,d", po::value<vector<int>>()->default_value({0}, "0"), "device ID of GPU (multiple times)")
			("chunkSize", po::value<int>()->default_value(1000), "number of concurrently processed structures at the server");
#endif
	options.push_back(concurrency);
	po::options_description sim("simulation");
	sim.add_options()
			("dielec", po::value<string>()->default_value("variable"), "dielectric behavior ('variable', 'constant')")
			("epsilon", po::value<double>()->default_value(15.0), "dielectric constant");
	options.push_back(sim);


	return options;
}

void TwoBodyParser::enforceRules(po::variables_map const& vm) const {
	std::vector<string> mutualExlusiveArgs = {"numCPUs", "device"};
	enforceMutualExcusiveness(vm, mutualExlusiveArgs);

	std::vector<string> allowedValues = {"single", "double"};
	enforceAllowedValues(vm, "prec", allowedValues);


	allowedValues = {"variable", "constant"};
	enforceAllowedValues(vm, "dielec", allowedValues);

	enforceUniqueness<int>(vm, "device");


}

void TwoBodyParser::assigneArgs(po::variables_map const& vm) noexcept {
	if(vm.count("dof"))
		_args->dofName = vm["dof"].as<string>();
	if(vm.count("receptor-pdb"))
		_args->recName = vm["receptor-pdb"].as<string>();
	if(vm.count("ligand-pdb"))
		_args->ligName = vm["ligand-pdb"].as<string>();
	if(vm.count("grid"))
		_args->gridName = vm["grid"].as<string>();
	if(vm.count("par"))
		_args->paramsName = vm["par"].as<string>();
	if(vm.count("alphabet"))
		_args->alphabetName = vm["alphabet"].as<string>();
	if(vm.count("numCPUs"))
		_args->numCPUs = vm["numCPUs"].as<int>();
	if(vm.count("device"))
		_args->deviceIds = vm["device"].as<vector<int>>();
	if(vm.count("chunkSize"))
		_args->chunkSize = vm["chunkSize"].as<int>();
	if(vm.count("prec"))
		_args->precision = vm["prec"].as<string>();
	if(vm.count("dielec"))
		_args->dielec = vm["dielec"].as<string>();
	if(vm.count("epsilon"))
		_args->epsilon = vm["epsilon"].as<double>();

}


