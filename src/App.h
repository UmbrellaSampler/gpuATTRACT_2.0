/*
 * App.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#ifndef SRC_APP_H_
#define SRC_APP_H_

#include "CmdArgs.h"

namespace as {

template<typename REAL>
class App {
	static_assert(std::is_arithmetic<REAL>::value, "Only arithmetic types supported");
public:
	virtual ~App() {};
	virtual void init(CmdArgs const& args) = 0;
	virtual void finalize() = 0;
	virtual void run() = 0;
};

}  // namespace as



#endif /* SRC_APP_H_ */
