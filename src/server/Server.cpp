/*
 * Server.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: uwe
 */

#include "ServerIncludes.h"
#include "CPU_6D_EnergyService.h"

#ifdef CUDA
#include "GPU_6D_EnergyService.h"
#endif

namespace as {
template
class Server<CPU_6D_EnergyService<float>>;

template
class Server<CPU_6D_EnergyService<double>>;

#ifdef CUDA
template
class Server<GPU_6D_EnergyService<float>>;

template
class Server<GPU_6D_EnergyService<double>>;
#endif

} // namespace as


