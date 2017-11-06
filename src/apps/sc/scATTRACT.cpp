/*
 * scATTRACT.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#include <CPUEnergyService6D.h>
#include <CPUEnergyService6DModes.h>
#include "scATTRACT.tpp"

#ifdef CUDA
#include <GPUEnergyService6D.h>
#endif

namespace as {

template
class scATTRACT<CPUEnergyService6D<float>>;

template
class scATTRACT<CPUEnergyService6D<double>>;

template
class scATTRACT<CPUEnergyService6DModes<float>>;

template
class scATTRACT<CPUEnergyService6DModes<double>>;

#ifdef CUDA
template
class scATTRACT<GPUEnergyService6D<float>>;

template
class scATTRACT<GPUEnergyService6D<double>>;
#endif

}  // namespace as

