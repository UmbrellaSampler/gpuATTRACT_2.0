/*
 * Server.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: uwe
 */

#include "ServerIncludes.h"
#include "Types_6D.h"
#include "Types_6D_Modes.h"
#include "Types_6D_MB_Modes.h"

namespace as {

template
class Server<Types_6D<float>>;

template
class Server<Types_6D<double>>;

template
class Server<Types_6D_Modes<float>>;

template
class Server<Types_6D_Modes<double>>;

template
class Server<Types_6D_MB_Modes<float>>;

template
class Server<Types_6D_MB_Modes<double>>;


} // namespace as


