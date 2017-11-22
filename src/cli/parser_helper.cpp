/*
 * CmdParserHelper.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include "parser_helper.h"
#include <sstream>
#include <numeric>

using namespace std;
namespace po = boost::program_options;

namespace as {
void enforceMutualExcusiveness(po::variables_map const& vm,
			vector<string> const& opts)
{
	for (size_t i = 0; i < opts.size()-1; ++i) {
		for (size_t j = i+1; j < opts.size(); j++) {
			if (vm.count(opts[i]) && !vm[opts[i]].defaulted() &&
				vm.count(opts[j]) && !vm[opts[j]].defaulted())
			{
				throw po::error(string("conflicting options '") +
									   opts[i] + "' and '" + opts[j] + "'.");
			}
		}
	}
}

template<typename T>
ostream& operator<< (ostream& s, const vector<T>& vec) {

	for (size_t i = 0; i < vec.size()-1; ++i ) {
		s << "'"<< vec[i] << "', ";
	}
	s << "'"<< *(vec.end()-1) << "'";
	return s;
}

template<typename T>
void enforceAllowedValues(po::variables_map const& vm,
		string opt, vector<T> const& values) {

	if (!vm.count(opt)) {
		throw logic_error("program option '" + opt + "' is not legal");
	}

	bool passed = false;
	for (auto const& val : values) {
		if (val == vm[opt].as<T>()) {
			passed = true;
		}
	}

	if (!passed) {
		stringstream msg;
		msg << "illegal option specification. option '" << opt << "' must be one of the following: " << values;
		throw po::error(msg.str());
	}
}

template<typename T>
void enforceUniqueness(po::variables_map const& vm, string opt)  {
	auto vec = vm[opt].as<vector<T>>();
	for (size_t i = 0; i < vec.size()-1; ++i) {
		for (size_t j = i+1; j < vec.size(); j++) {
			if (vec[i] == vec[j])
			{
				stringstream msg;
				msg << "illegal option specification. option '" << opt << "' is specified multiple times with the same value: " << vec[i];
				throw po::error(msg.str());
			}
		}
	}
}

string argumentUsage(po::options_description const& optsDesc) noexcept {
	(void) optsDesc;
//	const auto& opts = optsDesc.options();
//	stringstream usageMsg;
//	for (const auto& opt : opts) {
//		usageMsg << opt->canonical_display_name(0) << endl;
//		usageMsg << opt->canonical_display_name(1) << endl;
//		usageMsg << opt->canonical_display_name(2) << endl;
//		usageMsg << opt->long_name() << endl;
//		usageMsg << opt->format_name() << endl;
//		usageMsg << opt->format_parameter() << endl;
//		usageMsg << endl;
//	}
//	return usageMsg.str();

	// decided not to print anything but just options. can be changed here.
	return "[options]";
}

string extractFileName(string path) noexcept {
	const size_t last_slash_idx = path.find_last_of("\\/");
	if (string::npos != last_slash_idx)
	{
		path.erase(0, last_slash_idx + 1);
	}

	// Remove extension if present.
	const size_t period_idx = path.rfind('.');
	if (string::npos != period_idx)
	{
		path.erase(period_idx);
	}

	return path;
}

string concatenate(initializer_list<char const * const> const& values) noexcept {
	return std::accumulate(values.begin(), values.end(), std::string(),
    [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "', '" : "") + b;
    } );
}

std::string descriptionWithOptions(std::string desc, std::initializer_list<char const * const> const& values) noexcept {
	return string(desc)
			.append(" ('")
			.append(concatenate(values))
			.append("')");
}



template
void enforceAllowedValues(boost::program_options::variables_map const& vm,
		string opt, vector<string> const& values);

template
void enforceAllowedValues(boost::program_options::variables_map const& vm,
		string opt, vector<int> const& values);

template
void enforceUniqueness<int>(po::variables_map const& vm, string opt);

//template<typename T>
//void enforceUniqueness(po::variables_map const& vm, string opt);

} // namespace as
