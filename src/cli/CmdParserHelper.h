/*
 * CmdParserHelper.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_CMDPARSERHELPER_H_
#define SRC_CMDPARSERHELPER_H_

#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace as {

	void enforceMutualExcusiveness(boost::program_options::variables_map const& vm,
			std::vector<std::string> const& opts);

	template<typename T>
	void enforceAllowedValues(boost::program_options::variables_map const& vm,
			std::string opt, std::vector<T> const& values);

	template<typename T>
	void enforceUniqueness(boost::program_options::variables_map const& vm, std::string opt);
}



#endif /* SRC_CMDPARSERHELPER_H_ */
