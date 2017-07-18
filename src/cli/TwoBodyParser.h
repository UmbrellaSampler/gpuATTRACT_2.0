/*
 * TwoBodyParser.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_TWOBODYPARSER_H_
#define SRC_TWOBODYPARSER_H_

#include <boost/program_options.hpp>
#include <vector>

#include "AppCmdParser.h"

namespace as {

class TwoBodyParser : public AppCmdParser {
public:
	virtual ~TwoBodyParser() {};
	virtual void parse(int argc, char* argv[]) = 0;
protected:
	TwoBodyParser(std::shared_ptr<CmdArgs> args) : AppCmdParser::AppCmdParser(args) {}

	std::vector<boost::program_options::options_description> options() const noexcept;
	void enforceRules(boost::program_options::variables_map const& vm) const;
	void assigneArgs(boost::program_options::variables_map const& vm) noexcept;

};

}



#endif /* SRC_TWOBODYPARSER_H_ */
