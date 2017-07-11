/*
 * RequestHandler.cpp
 *
 *  Created on: Jul 9, 2017
 *      Author: uwe
 */


#include "RequestHandler.tpp"
#include "Server.h"
#include "Request.h"
#include "CPU_6D_EnergyService.h"

#ifdef CUDA
#include "GPU_6D_EnergyService.h"
#endif

namespace as {

template
class RequestHandler<Server<CPU_6D_EnergyService<float>>>;

template
class RequestHandler<Server<CPU_6D_EnergyService<double>>>;

#ifdef CUDA
template
class RequestHandler<Server<GPU_6D_EnergyService<float>>>;

template
class RequestHandler<Server<GPU_6D_EnergyService<double>>>;
#endif

}  // namespace as

