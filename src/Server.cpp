/*
 * Server.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: uwe
 */

#include "ServerIncludes.h"
#include "CPU_6D_EnergyService.h"

namespace as {
template
class Server<CPU_6D_EnergyService<float>>;

template
class Server<CPU_6D_EnergyService<double>>;
}


