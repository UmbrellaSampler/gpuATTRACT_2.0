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

protected:
	TwoBodyParser(std::shared_ptr<CmdArgs> args) : AppCmdParser::AppCmdParser(args) {}

	virtual void addOptions() noexcept override;
	void enforceRules(boost::program_options::variables_map const& vm) const override;
	void assignArgs(boost::program_options::variables_map const& vm) noexcept override;

};

}



#endif /* SRC_TWOBODYPARSER_H_ */
