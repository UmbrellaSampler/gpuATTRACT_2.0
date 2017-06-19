/*
 * ScoreParser.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_SCOREPARSER_H_
#define SRC_SCOREPARSER_H_

#include "TwoBodyParser.h"

namespace as {

class ScoreParser : public TwoBodyParser {
public:
	ScoreParser(std::shared_ptr<CmdArgs> args) : TwoBodyParser::TwoBodyParser(args) {}

	std::vector<boost::program_options::options_description> options();
	void parse(int argc, char* argv[]) override;
	void enforceRules(boost::program_options::variables_map const& vm) const;
	void assigneArgs(boost::program_options::variables_map const& vm) noexcept;
};

}  // namespace as



#endif /* SRC_SCOREPARSER_H_ */
