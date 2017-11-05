/*
 * AppCmdParser.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: uwe
 */

#include <parser_helper.h>
#include <iostream>
#include <fstream>
#include "AppCmdParser.h"

namespace po = boost::program_options;

namespace as {

void AppCmdParser::parse(int argc, char* argv[])  {
	_applicationName = extractFileName(argv[0]);
	addOptions();
	po::variables_map vm;
	try {
		store(po::command_line_parser(std::min(argc, argc), argv).
				  options(_optsDesc).run(), vm);

		if (vm.count("help")) {
			printHelp();
			exit(EXIT_SUCCESS);
		}

		notify(vm);

		if (vm.count("config")) {
			std::string config_file = vm["config"].as<std::string>();
			std::ifstream ifs(config_file.c_str());
			if (!ifs)
			{
				throw po::error("cannot open config file: " + config_file + "\n");
			}
			else
			{
				store(parse_config_file(ifs, _optsDesc), vm);
				notify(vm);
			}
		}
		enforceRules(vm);
		assignArgs(vm);

	} catch (po::error& e) {
		std::cerr << "error: " << e.what() << std::endl;
		printHelp();
		exit(EXIT_FAILURE);
	} catch (std::exception& e) {
		std::cerr << "unexpected exception after cmd-line parsing: " << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
}

void AppCmdParser::printHelp() const noexcept {
	std::cout << std::endl << usage() << std::endl;
	std::cout << _optsDesc << std::endl;
}

std::string AppCmdParser::usage() const noexcept {
	return "usage: " + _applicationName + " " + argumentUsage(_optsDesc);
}

}
