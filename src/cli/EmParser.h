/*
 * EmParser.h
 *
 *  Created on: Nov 4, 2017
 *      Author: uwe
 */

#ifndef SRC_CLI_EMPARSER_H_
#define SRC_CLI_EMPARSER_H_

#include "TwoBodyParser.h"

namespace as {

class EmParser : public TwoBodyParser {
public:
	EmParser(std::shared_ptr<CmdArgs> args) : TwoBodyParser::TwoBodyParser(args) {}
	virtual ~EmParser() {};

protected:
	std::string appShortName() const noexcept override;
	void addOptions() noexcept override;
	void enforceRules(boost::program_options::variables_map const& vm) const override;
	void assignArgs(boost::program_options::variables_map const& vm) noexcept override;
};

}  // namespace as



#endif /* SRC_CLI_EMPARSER_H_ */
