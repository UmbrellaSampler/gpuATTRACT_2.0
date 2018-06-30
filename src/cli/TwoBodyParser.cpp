/*
 * TwoBodyParser.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */


#include "parser_helper.h"
#include "TwoBodyParser.h"
#include "parser_constants.h"

#include <string>

using namespace as;
using namespace std;
namespace po = boost::program_options;

void TwoBodyParser::addOptions() noexcept {

	po::options_description config("generic");
	config.add_options()
			("help", "print this help message")
			("config", po::value<string>(), "configuration file")
			("prec", po::value<string>()->default_value("single"),
					descriptionWithOptions("arithmetic precision", GENERIC_ALLOWED_PRECISION).c_str());
	_optsDesc.add(config);

	po::options_description input("input files");
	input.add_options()
			("dof"     			  , po::value<string>()->required()									, "structure (DOF) file")
			("receptor-pdb,r"     , po::value<string>()->default_value(FILE_DEFAULT_RECEPTOR_PDB)	, "pdb-file of receptor")
			("ligand-pdb,l"       , po::value<string>()->default_value(FILE_DEFAULT_LIGANG_PDB)   	, "pdb-file of ligand")
			("gridrec"             , po::value<string>()->default_value(FILE_DEFAULT_RECEPTOR_GRID)	, "receptor grid file")
			("gridlig"             , po::value<string>()->default_value(FILE_DEFAULT_LIGAND_GRID)	, "ligand grid file")
			("par,p"	          , po::value<string>()->default_value(FILE_DEFAULT_PARAMETER)		, "attract forcefield parameter file")
			("alphabetrec"		  , po::value<string>()->default_value(FILE_DEFAULT_GRID_ALPAHBET_RECEPTOR)	, "receptor grid alphabet file")
			("alphabetlig"		  , po::value<string>()->default_value(FILE_DEFAULT_GRID_ALPAHBET_LIGAND)	, "ligand grid alphabet file")
			("modesl"	          , po::value<string>()->default_value(DEFAULT_MODE_RECEPTOR_FILE)  , "mode file of ligand")
			("modesr"	          , po::value<string>()->default_value(DEFAULT_MODE_LIGAND_FILE)	, "mode file of receptor");
	_optsDesc.add(input);

	po::options_description concurrency("concurrency");

	concurrency.add_options()
#ifndef CUDA
			("numCPUs,c", po::value<int>()->default_value(SERVER_DEFAULT_NUM_CPUS), "number of CPU threads for CPU mode")
#else
			("numCPUs,c", po::value<int>()->default_value(0), "number of CPU threads for CPU mode")
			("device,d", po::value<vector<int>>()->default_value({SERVER_DEFAULT_DEVICE_ID}, std::to_string(SERVER_DEFAULT_DEVICE_ID)), "device ID of GPU (multiple times)")
#endif
			("chunkSize", po::value<int>()->default_value(SERVER_DEFAULT_CHUNK_SIZE), "number of concurrently processed structures at the server");

	_optsDesc.add(concurrency);
	po::options_description sim("simulation");
	sim.add_options()
		("numModes", po::value<int>()->default_value(DEFAULT_NUM_MODES), "number of modes")
		("dielec", po::value<string>()->default_value(SIM_DEFAULT_DIELEC),	descriptionWithOptions("dielectric behavior", SIM_ALLOWED_DIELEC).c_str())
		("epsilon", po::value<double>()->default_value(SIM_DEFAULT_EPSILON), "dielectric constant")
		("evscale", po::value<double>()->default_value(SIM_DEFAULT_EVSCALE), "scales the eigenvalues of the modes. Default is 1.0.")
		("cutoff", po::value<double>()->default_value(SIM_DEFAULT_CUTOFF), "squared cutoff radius. default switched off( >10000)");
	_optsDesc.add(sim);

}

void TwoBodyParser::enforceRules(po::variables_map const& vm) const {
	std::vector<string> mutualExlusiveArgs = {"numCPUs", "device"};
	enforceMutualExcusiveness(vm, mutualExlusiveArgs);

	enforceAllowedValues(vm, "prec", vector<string>(GENERIC_ALLOWED_PRECISION.begin(),
					GENERIC_ALLOWED_PRECISION.end()));

	enforceAllowedValues(vm, "dielec", vector<string>(SIM_ALLOWED_DIELEC.begin(),
					SIM_ALLOWED_DIELEC.end()));

#ifdef CUDA
	enforceUniqueness<int>(vm, "device");
#endif


}

void TwoBodyParser::assignArgs(po::variables_map const& vm) noexcept {
	if(vm.count("dof"))
		_args->dofName = vm["dof"].as<string>();
	if(vm.count("receptor-pdb"))
		_args->recName = vm["receptor-pdb"].as<string>();
	if(vm.count("ligand-pdb"))
		_args->ligName = vm["ligand-pdb"].as<string>();
	if(vm.count("gridrec"))
		_args->gridRecName = vm["gridrec"].as<string>();
	if(vm.count("gridlig"))
		_args->gridLigName = vm["gridlig"].as<string>();
	if(vm.count("par"))
		_args->paramsName = vm["par"].as<string>();
	if(vm.count("alphabetrec"))
		_args->alphabetRecName = vm["alphabetrec"].as<string>();
	if(vm.count("alphabetlig"))
		_args->alphabetLigName = vm["alphabetlig"].as<string>();
	if(vm.count("modesr"))
		_args->recModesName = vm["modesr"].as<string>();
	if(vm.count("modesl"))
		_args->ligModesName = vm["modesl"].as<string>();
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

	if(vm.count("cutoff"))
			_args->cutoff = vm["cutoff"].as<double>();
	if(vm.count("evscale"))
			_args->evscale = vm["evscale"].as<double>();
	if(vm.count("numModes"))
		_args->numModes = vm["numModes"].as<int>();

}


