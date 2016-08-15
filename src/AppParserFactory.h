/*
 * AppParserFactory.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_APPPARSERFACTORY_H_
#define SRC_APPPARSERFACTORY_H_

#include <memory>
#include "AppCmdParser.h"

namespace as {

class AppParserFactory {
public:
	static std::unique_ptr<AppCmdParser> create(AppType app, std::shared_ptr<CmdArgs> args);
};

}  // namespace as

#endif /* SRC_APPPARSERFACTORY_H_ */
