/*
 * Server.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: uwe
 */

#include "ServerIncludes.h"
#include "Types_6D.h"

namespace as {

template
class Server<Types_6D<float>>;

template
class Server<Types_6D<double>>;


} // namespace as


