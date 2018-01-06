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
			("dofRec"     			  , po::value<string>()->required()									, "receptor structure (DOF) file")
			("dofLig"     			  , po::value<string>()->required()									, "ligand structure (DOF) file")
			("dofLig2"     			  , po::value<string>()->required()									, "second ligand structure (DOF) file")
			("dofLig3"     			  , po::value<string>()->required()									, "third ligand structure (DOF) file")
			("receptor-pdb,r"     , po::value<string>()->default_value(FILE_DEFAULT_RECEPTOR_PDB)	, "pdb-file of receptor")
			("ligand-pdb,l"       , po::value<string>()->default_value(FILE_DEFAULT_LIGANG_PDB)   	, "pdb-file of ligand")
			("ligand-pdb2,l"       , po::value<string>()->default_value(FILE_DEFAULT_LIGANG_PDB)   	, "second pdb-file of ligand")
			("ligand-pdb3,l"       , po::value<string>()->default_value(FILE_DEFAULT_LIGANG_PDB)   	, "third pdb-file of ligand")
			("gridrec,gr"             , po::value<string>()->default_value(FILE_DEFAULT_RECEPTOR_GRID)	, "receptor grid file")
			("gridlig,gl"             , po::value<string>()->default_value(FILE_DEFAULT_LIGAND_GRID)	, "ligand grid file")
			("gridlig2,gl2"             , po::value<string>()->default_value(FILE_DEFAULT_LIGAND_GRID)	, "second ligand grid file")
			("gridlig3,gl3"             , po::value<string>()->default_value(FILE_DEFAULT_LIGAND_GRID)	, "third ligand grid file")
			("par,p"	          , po::value<string>()->default_value(FILE_DEFAULT_PARAMETER)		, "attract forcefield parameter file")
			("alphabetrec,ar"		  , po::value<string>()->default_value(FILE_DEFAULT_GRID_ALPAHBET_RECEPTOR)	, "receptor grid alphabet file")
			("alphabetlig,al"		  , po::value<string>()->default_value(FILE_DEFAULT_GRID_ALPAHBET_LIGAND)	, "ligand grid alphabet file")
			("modl,ml"	          , po::value<string>()->default_value(DEFAULT_MODE_RECEPTOR_FILE)  , "mode file of ligand")
			("modl2,ml"	          , po::value<string>()->default_value(DEFAULT_MODE_RECEPTOR_FILE)  , "second mode file of ligand")
			("modl3,ml"	          , po::value<string>()->default_value(DEFAULT_MODE_RECEPTOR_FILE)  , "third mode file of ligand")
			("modr,mr"	          , po::value<string>()->default_value(DEFAULT_MODE_LIGAND_FILE)	, "mode file of receptor")
			("numLig,mb"	          , po::value<string>()->default_value(DEFAULT_NUM_LIGANDS)	, "number of Ligands");
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
			("dielec", po::value<string>()->default_value(SIM_DEFAULT_DIELEC),
					descriptionWithOptions("dielectric behavior", SIM_ALLOWED_DIELEC).c_str())
			("epsilon", po::value<double>()->default_value(SIM_DEFAULT_EPSILON), "dielectric constant");
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
	if(vm.count("dofRec"))
		_args->dofName = vm["dofRec"].as<string>();
	if(vm.count("dofLig"))
		_args->dofNameLig.push_back(vm["dofLig"].as<string>());
	if(vm.count("dofLig2"))
		_args->dofNameLig.push_back(vm["dofLig2"].as<string>());
	if(vm.count("dofLig3"))
		_args->dofNameLig.push_back(vm["dofLig3"].as<string>());
	if(vm.count("receptor-pdb"))
		_args->recName = vm["receptor-pdb"].as<string>();
	if(vm.count("ligand-pdb"))
		_args->ligName.push_back(vm["ligand-pdb"].as<string>());
	if(vm.count("ligand-pdb2"))
		_args->ligName.push_back(vm["ligand-pdb2"].as<string>());
	if(vm.count("ligand-pdb3"))
		_args->ligName.push_back(vm["ligand-pdb3"].as<string>());
	if(vm.count("gridrec"))
		_args->gridRecName = vm["gridrec"].as<string>();
	if(vm.count("gridlig"))
		_args->gridLigName.push_back(vm["gridlig"].as<string>());
	if(vm.count("gridlig2"))
		_args->gridLigName.push_back(vm["gridlig2"].as<string>());
	if(vm.count("gridlig3"))
		_args->gridLigName.push_back(vm["gridlig3"].as<string>());
	if(vm.count("par"))
		_args->paramsName = vm["par"].as<string>();
	if(vm.count("alphabetrec"))
		_args->alphabetRecName = vm["alphabetrec"].as<string>();
	if(vm.count("alphabetlig"))
		_args->alphabetLigName = vm["alphabetlig"].as<string>();
	if(vm.count("modr"))
		_args->recModesName = vm["modr"].as<string>();
	if(vm.count("modl"))
		_args->ligModesName.push_back(vm["modl"].as<string>());
	if(vm.count("modl2"))
		_args->ligModesName.push_back(vm["modl2"].as<string>());
	if(vm.count("modl3"))
		_args->ligModesName.push_back(vm["modl3"].as<string>());
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
	if(vm.count("numLigands"))
		_args->numLigands = vm["numLigands"].as<int>();

}


