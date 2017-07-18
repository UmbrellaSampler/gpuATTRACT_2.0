/*
 * TwoBodyParser.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <iostream>
#include <string>
#include <fstream>
#include "TwoBodyParser.h"
#include "CmdParserHelper.h"

using namespace as;
using namespace std;
namespace po = boost::program_options;

void TwoBodyParser::parse(int argc, char* argv[])  {
	std::vector<po::options_description> opts = options();
	po::variables_map vm;
	try {
		po::options_description cmdline_options;
		for (auto const& opt: opts) {
			cmdline_options.add(opt);
		}

		store(po::command_line_parser(std::min(argc, argc), argv).
				  options(cmdline_options).run(), vm);


		if (vm.count("help")) {
			cout << usage() << endl << endl;
			cout << cmdline_options << endl;
			exit(EXIT_SUCCESS);
		}

		notify(vm);

		if (vm.count("config")) {
			string config_file = vm["config"].as<string>();
			ifstream ifs(config_file.c_str());
			if (!ifs)
			{
				throw po::error("cannot open config file: " + config_file + "\n");
			}
			else
			{
				store(parse_config_file(ifs, cmdline_options), vm);
				notify(vm);
			}
		}
		enforceRules(vm);
		assigneArgs(vm);

	} catch (po::error& e) {
		cerr << "error: " << e.what() << endl;
		cout << usage() << "\n\n";
		exit(EXIT_FAILURE);
	} catch (std::exception& e) {
		cerr << "unexpected exception after cmd-line parsing: " << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

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
			("ligand-pdb,l"       , po::value<string>()->default_value("ligandr.pdb")   		, "pdb-file of ligand")
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

string TwoBodyParser::usage() const noexcept {
	stringstream msg;
	msg << "usage: gpuAttract sc --dof <file> [--config <file>] [-r <file>] [-l <file>] [-g <file>] [-p <file>] [-a <file>] [--prec <string>]";
	msg << "[-c <int>] ";
#ifdef CUDA
	msg << "[-d <int>...] ";
#endif
	msg << "[--chunkSize <int>]";
	return msg.str();
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


