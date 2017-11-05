/*
 * AppFactory.h
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#ifndef SRC_APPFACTORY_H_
#define SRC_APPFACTORY_H_


#include <memory>
#include "AppType.h"

namespace as {

class App;
class CmdArgs;

class AppFactory {
public:
	AppFactory() = delete;
	static std::unique_ptr<App> create(const CmdArgs& args);

private:
	template<typename REAL>
	static std::unique_ptr<App> create(AppType appType);
};

}  // namespace as


#endif /* SRC_APPFACTORY_H_ */
