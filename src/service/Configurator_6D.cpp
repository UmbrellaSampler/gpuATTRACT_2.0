/*
 * Init_6D.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#include "Configurator_6D.tpp"
//#include "CPU_6D_EnergyService.h"
//#ifdef CUDA
//#include "GPU_6D_EnergyService.h"
//#endif

namespace as {

template
class Configurator_6D<float>;

template
class Configurator_6D<double>;

}  // namespace as

