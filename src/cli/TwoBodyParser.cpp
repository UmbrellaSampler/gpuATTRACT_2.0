/*
 * TwoBodyParser.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */


#include <parser_helper.h>
#include "TwoBodyParser.h"
#include "parser_constants.h"

using namespace as;
using namespace std;
namespace po = boost::program_options;

void TwoBodyParser::addOptions() noexcept {

	po::options_description config("generic");
	config.add_options()
			("help", "print this help message")
			("config", po::value<string>(), "configuration file")
			("prec", po::value<string>()->default_value("single"), "arithmetic precision ('single', 'double')");
	_optsDesc.add(config);

	po::options_description input("input files");
	input.add_options()
			("dof"     			  , po::value<string>()->required()									, "structure (DOF) file")
			("receptor-pdb,r"     , po::value<string>()->default_value(DEFAULT_RECEPTOR_PDB_FILE)	, "pdb-file of receptor")
			("ligand-pdb,l"       , po::value<string>()->default_value(DEFAULT_LIGANG_PDB_FILE)   	, "pdb-file of ligand")
			("grid,g"             , po::value<string>()->default_value(DEFAULT_RECEPTOR_GRID_FILE)	, "receptor grid file")
			("par,p"	          , po::value<string>()->default_value(DEFAULT_PARAMETER_FILE)		, "attract forcefield parameter file")
			("alphabet,a"		  , po::value<string>()->default_value(DEFAULT_GRID_ALPAHBET_FILE)	, "receptor grid alphabet file");
			("modl,ml"	          , po::value<string>()->default_value(DEFAULT_MODE_RECEPTOR_FILE)  , "mode file of ligand");
			("modr,mr"	          , po::value<string>()->default_value(DEFAULT_MODE_LIGAND_FILE)	, "mode file of receptor");
	_optsDesc.add(input);

	po::options_description concurrency("concurrency");

	concurrency.add_options()
#ifndef CUDA
			("numCPUs,c", po::value<int>()->default_value(DEFAULT_NUM_CPUS), "number of CPU threads for CPU mode")
#else
			("numCPUs,c", po::value<int>()->default_value(0), "number of CPU threads for CPU mode")
			("device,d", po::value<vector<int>>()->default_value({0}, "0"), "device ID of GPU (multiple times)")
#endif
			("chunkSize", po::value<int>()->default_value(DEFAULT_CHUNK_SIZE), "number of concurrently processed structures at the server");

	_optsDesc.add(concurrency);
	po::options_description sim("simulation");
	sim.add_options()
			("numModes", po::value<int>()->default_value(DEFAULT_NUM_MODES), "number of modes")
			("dielec", po::value<string>()->default_value(DEFAULT_DIELEC_MODE), "dielectric behavior ('variable', 'constant')")
			("epsilon", po::value<double>()->default_value(DEFAULT_EPSILON_CONSTANT), "dielectric constant");
	_optsDesc.add(sim);

}

void TwoBodyParser::enforceRules(po::variables_map const& vm) const {
	std::vector<string> mutualExlusiveArgs = {"numCPUs", "device"};
	enforceMutualExcusiveness(vm, mutualExlusiveArgs);

	std::vector<string> allowedValues = {"single", "double"};
	enforceAllowedValues(vm, "prec", allowedValues);


	allowedValues = {"variable", "constant"};
	enforceAllowedValues(vm, "dielec", allowedValues);

#ifdef CUDA
	enforceUniqueness<int>(vm, "device");
#endif


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
	if(vm.count("modr"))
				_args->recModesName = vm["modr"].as<string>();
	if(vm.count("modl"))
			_args->ligModesName = vm["modl"].as<string>();
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
	if(vm.count("numModes"))
		_args->numModes = vm["numModes"].as<int>();

}


