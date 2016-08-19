/*
 * CmdParser.h
 *
 *  Created on: Aug 14, 2016
 *      Author: uwe
 */

#ifndef SRC_CMDPARSER_H_
#define SRC_CMDPARSER_H_

#include <string>
#include <memory>

#include "CmdArgs.h"

namespace as {

class CmdParser {
public:
	CmdParser();

	void parse(int argc, char* argv[]) noexcept;

	CmdArgs args() {
		return *_args;
	}
private:

	void assignApp(std::string);

	std::shared_ptr<CmdArgs> _args;
};

}  // namespace as



#endif /* SRC_CMDPARSER_H_ */
