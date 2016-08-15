/*
 * CmdParserHelper.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <sstream>
#include "CmdParserHelper.h"

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
void enforceAllowedValues(boost::program_options::variables_map const& vm,
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
		msg << "illegal option assignment. option '" << opt << "' must be one of the following: " << values;
		throw po::error(msg.str());
	}
}

template
void enforceAllowedValues(boost::program_options::variables_map const& vm,
		std::string opt, std::vector<std::string> const& values);

template
void enforceAllowedValues(boost::program_options::variables_map const& vm,
		std::string opt, std::vector<int> const& values);

} // namespace as
