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

	std::string argumentUsage(boost::program_options::options_description const& optsDesc) noexcept;

	std::string extractFileName(std::string path) noexcept;

	std::string concatenate(std::initializer_list<char const * const> const& values) noexcept;
	std::string descriptionWithOptions(std::string desc, std::initializer_list<char const * const> const& values) noexcept;

}



#endif /* SRC_CMDPARSERHELPER_H_ */
