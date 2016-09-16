/*
 * Init_6D.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#include "Configurator_6D.tpp"
#include "CPU_6D_EnergyService.h"
#ifdef CUDA
#include "GPU_6D_EnergyService.h"
#endif

namespace as {

template
class Configurator_6D<CPU_6D_EnergyService<float>>;

template
class Configurator_6D<CPU_6D_EnergyService<double>>;

#ifdef CUDA
template
class Configurator_6D<GPU_6D_EnergyService<float>>;

template
class Configurator_6D<GPU_6D_EnergyService<double>>;
#endif

}  // namespace as

