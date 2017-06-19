/*
 * AppCmdParser.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_APPCMDPARSER_H_
#define SRC_APPCMDPARSER_H_

#include <memory>
#include "CmdArgs.h"

namespace as {

class AppCmdParser {
public:
	AppCmdParser(std::shared_ptr<CmdArgs> args) : _args(args) {}
	virtual ~AppCmdParser() {};

	virtual void parse(int argc, char* argv[]) = 0;

protected:
	std::shared_ptr<CmdArgs> _args;
};

} // namespace



#endif /* SRC_APPCMDPARSER_H_ */
