/*
 * ScoreParser.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_MCPARSER_H_
#define SRC_MCPARSER_H_

#include "TwoBodyParser.h"

namespace as {

class McParser : public TwoBodyParser {
public:
	McParser(std::shared_ptr<CmdArgs> args) : TwoBodyParser::TwoBodyParser(args) {}
	virtual ~McParser() {};

protected:
	std::string appShortName() const noexcept override;

};

}  // namespace as



#endif /* SRC_SCOREPARSER_H_ */
