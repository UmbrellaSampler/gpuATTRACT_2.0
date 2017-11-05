/*
 * EmParser.cpp
 *
 *  Created on: Nov 4, 2017
 *      Author: uwe
 */

#include "parser_helper.h"
#include "EmParser.h"
#include "parser_constants.h"

using namespace std;
namespace po = boost::program_options;

namespace as {

std::string EmParser::appShortName() const noexcept {
	return APP_SHORT_NAME_EM;
}

void EmParser::addOptions() noexcept {
	TwoBodyParser::addOptions();
	po::options_description minimization("minimization");
	minimization.add_options()
			("solver,s", po::value<string>()->default_value(EM_DEFAULT_SOLVER),
			descriptionWithOptions("Optimization algorithms", EM_ALLOWED_SOLVERS).c_str())
			("maxConcurrency", po::value<int>()->default_value(EM_DEFAULT_CONCURRENCY),
			"Max. number of concurrent structures that may be processed at the same time")
			("numChunks", po::value<int>()->default_value(EM_DEFAULT_NUM_CHUNKS),
			"Number of request chunks");
	_optsDesc.add(minimization);
}

void EmParser::enforceRules(
		boost::program_options::variables_map const& vm) const {
	TwoBodyParser::enforceRules(vm);
	enforceAllowedValues(vm, "solver",
			vector<string>(EM_ALLOWED_SOLVERS.begin(),
					EM_ALLOWED_SOLVERS.end()));
}

void EmParser::assignArgs(po::variables_map const& vm) noexcept {
	TwoBodyParser::assignArgs(vm);
	if(vm.count("solver"))
		_args->solver = vm["solver"].as<string>();
	if(vm.count("maxConcurrency"))
		_args->maxConcurrency = vm["maxConcurrency"].as<int>();
	if(vm.count("numChunks"))
		_args->numChunks = vm["numChunks"].as<int>();
}

}

