/*
 * CmdParserHelper.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <parser_helper.h>
#include <sstream>

using namespace std;
namespace po = boost::program_options;

namespace as {
void enforceMutualExcusiveness(po::variables_map const& vm,
			std::vector<std::string> const& opts)
{
	for (size_t i = 0; i < opts.size()-1; ++i) {
		for (size_t j = i+1; j < opts.size(); j++) {
			if (vm.count(opts[i]) && !vm[opts[i]].defaulted() &&
				vm.count(opts[j]) && !vm[opts[j]].defaulted())
			{
				throw po::error(std::string("conflicting options '") +
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
		std::string opt, std::vector<T> const& values) {

	if (!vm.count(opt)) {
		throw std::logic_error("program option '" + opt + "' is not legal");
	}

	bool passed = false;
	for (auto const& val : values) {
		if (val == vm[opt].as<T>()) {
			passed = true;
		}
	}

	if (!passed) {
		std::stringstream msg;
		msg << "illegal option specification. option '" << opt << "' must be one of the following: " << values;
		throw po::error(msg.str());
	}
}

template<typename T>
void enforceUniqueness(po::variables_map const& vm, std::string opt)  {
	auto vec = vm[opt].as<vector<T>>();
	for (size_t i = 0; i < vec.size()-1; ++i) {
		for (size_t j = i+1; j < vec.size(); j++) {
			if (vec[i] == vec[j])
			{
				std::stringstream msg;
				msg << "illegal option specification. option '" << opt << "' is specified multiple times with the same value: " << vec[i];
				throw po::error(msg.str());
			}
		}
	}
}

std::string argumentUsage(po::options_description const& optsDesc) noexcept {
	(void) optsDesc;
//	const auto& opts = optsDesc.options();
//	std::stringstream usageMsg;
//	for (const auto& opt : opts) {
//		usageMsg << opt->canonical_display_name(0) << std::endl;
//		usageMsg << opt->canonical_display_name(1) << std::endl;
//		usageMsg << opt->canonical_display_name(2) << std::endl;
//		usageMsg << opt->long_name() << std::endl;
//		usageMsg << opt->format_name() << std::endl;
//		usageMsg << opt->format_parameter() << std::endl;
//		usageMsg << std::endl;
//	}
//	return usageMsg.str();

	// decided not to print anything but just options. can be changed here.
	return "[options]";
}

std::string extractFileName(std::string path) noexcept {
	const size_t last_slash_idx = path.find_last_of("\\/");
	if (std::string::npos != last_slash_idx)
	{
		path.erase(0, last_slash_idx + 1);
	}

	// Remove extension if present.
	const size_t period_idx = path.rfind('.');
	if (std::string::npos != period_idx)
	{
		path.erase(period_idx);
	}

	return path;
}



template
void enforceAllowedValues(boost::program_options::variables_map const& vm,
		std::string opt, std::vector<std::string> const& values);

template
void enforceAllowedValues(boost::program_options::variables_map const& vm,
		std::string opt, std::vector<int> const& values);

template
void enforceUniqueness<int>(po::variables_map const& vm, std::string opt);

//template<typename T>
//void enforceUniqueness(po::variables_map const& vm, std::string opt);

} // namespace as
