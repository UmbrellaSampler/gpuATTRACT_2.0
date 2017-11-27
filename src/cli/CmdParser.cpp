/*
 * CmdParser.cpp
 *
 *  Created on: Aug 14, 2016
 *      Author: uwe
 */

#include <boost/program_options.hpp>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "CmdParser.h"
#include "ScoreParser.h"
#include "version.h"
#include "AppParserFactory.h"
#include "parser_constants.h"

using namespace as;
using namespace std;
namespace po = boost::program_options;

CmdParser::CmdParser() {
	_args = std::make_shared<CmdArgs>();
}

static po::options_description appOptions () noexcept {
	po::options_description options("generic");
	options.add_options()
			("version,v", "print version information")
			("help", "print this help massage");

	return options;

}

static po::positional_options_description posOptions() {
	po::positional_options_description options;
	options.add("app", 1);
	return options;
}

static po::options_description hiddenOptions() {
	po::options_description options("Hidden");
	options.add_options()
			("app", po::value<string>(), "application: sc,scmode,em,mc");
	return options;
}

static string usage() noexcept {
	stringstream msg;
	msg << "usage: gpuAttract <app> <args>" << endl;
	msg << endl;
	msg << "available apps:" << endl;
	msg << "\t" << APP_SHORT_NAME_SC << "\t" << "Scoring input structures" << endl;
	msg << endl;
	msg << "\t" << APP_SHORT_NAME_EM << "\t" << "Docking of input structures using energy minimization" << endl;
	msg << endl;
	msg << "For further help information type: \"gpuAttract --help\" or \"gpuAttract <app> --help\"";
	return msg.str();
}

void CmdParser::parse(int argc, char* argv[]) noexcept {

	// get only the first argument to determine the application
	auto generic = appOptions();
	auto positional = posOptions();
	auto hidden = hiddenOptions();

	po::variables_map vm;
	try {
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(hidden);

		po::options_description visible;
		visible.add(generic);

		store(po::command_line_parser(std::min(2,argc), argv).
				  options(cmdline_options).positional(positional).run(), vm);

		if (vm.count("help")) {
			cout << usage() << endl << endl;
			cout << generic << endl;
			exit(EXIT_SUCCESS);
		}

		notify(vm);

		if (vm.empty()) {
			throw po::error("no argument specified");
		}

		if (vm.count("version")) {
			cout << "gpuATTRACT version " << MAJOR_VER << "." << MINOR_VER << endl;
		}

		if (vm.count("app")) {
			string app = vm["app"].as<string>();
			assignApp(app);
			std::unique_ptr<AppCmdParser> appParser = AppParserFactory::create(_args->app, _args);
			int argc_app = argc-1;
			char* argv_app[argc_app];
			argv_app[0] = argv[0];
			std::copy(argv+2, argv+argc, argv_app+1);
			appParser->parse(argc_app, argv_app);
		}

	} catch (po::error& e) {
		cerr << "error: " << e.what() << endl;
		cout << usage() << "\n\n";
		exit(EXIT_FAILURE);
	} catch (std::exception& e) {
		cerr << "unexpected exception after cmd-line parsing: " << e.what() << endl;
	}

}

void CmdParser::assignApp(string app) {
	if (app.compare(APP_SHORT_NAME_SC) == 0	) {
		_args->app = AppType::SCORE;
	} else if (app.compare(APP_SHORT_NAME_EM) == 0	) {
		_args->app = AppType::EM;
	} else if (app.compare(APP_SHORT_NAME_SCMODE) == 0	) {
		_args->app = AppType::SCORE_MODE;
	} else {
		throw po::error("unknown app: " + app);
	}
}

