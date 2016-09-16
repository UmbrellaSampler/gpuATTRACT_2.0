/*
 * scATTRACT.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#include "CPU_6D_EnergyService.h"
#include "scATTRACT.tpp"

#ifdef CUDA
#include "GPU_6D_EnergyService.h"
#endif

namespace as {

template
class scATTRACT<CPU_6D_EnergyService<float>>;

template
class scATTRACT<CPU_6D_EnergyService<double>>;

#ifdef CUDA
template
class scATTRACT<GPU_6D_EnergyService<float>>;

template
class scATTRACT<GPU_6D_EnergyService<double>>;
#endif

}  // namespace as

