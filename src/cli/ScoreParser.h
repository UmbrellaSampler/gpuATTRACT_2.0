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
	virtual ~ScoreParser() {};

protected:

};

}  // namespace as



#endif /* SRC_SCOREPARSER_H_ */
