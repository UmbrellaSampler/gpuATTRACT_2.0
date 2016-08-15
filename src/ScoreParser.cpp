/*
 * ScoreParser.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <iostream>
#include <sstream>
#include <fstream>

#include "ScoreParser.h"

using namespace as;
using namespace std;
namespace po = boost::program_options;

std::vector<po::options_description> ScoreParser::options() {
	return TwoBodyParser::options();
}

static string usage() noexcept {
	stringstream msg;
	msg << "usage: gpuAttract sc --dof <file> [--config <file>] [-r <file>] [-l <file>] [-g <file>] [-p <file>] [-a <file>] [--prec <string>]";
	msg << "[-c <int>] ";
#ifdef CUDA
	msg << "[-d <int>...] ";
#endif
	msg << "[--chunkSize <int>]";
	return msg.str();
}

void ScoreParser::parse(int argc, char* argv[])  {
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

void ScoreParser::enforceRules(po::variables_map const& vm) const {
	TwoBodyParser::enforceRules(vm);
}

void ScoreParser::assigneArgs(po::variables_map const& vm) noexcept {
	TwoBodyParser::assigneArgs(vm);
}


