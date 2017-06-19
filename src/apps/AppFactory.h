/*
 * AppFactory.h
 *
 *  Created on: Sep 4, 2016
 *      Author: uwe
 */

#ifndef SRC_APPFACTORY_H_
#define SRC_APPFACTORY_H_


#include <memory>
#include "CmdArgs.h"

namespace as {

class App;

class AppFactory {
public:
	static std::unique_ptr<App> create(const CmdArgs& args);

private:
	enum class Platform {
		CPU,
		GPU,
		unspecified
	};

	template<typename REAL>
	static std::unique_ptr<App> create(AppType appType, Platform p);



};

}  // namespace as


#endif /* SRC_APPFACTORY_H_ */
