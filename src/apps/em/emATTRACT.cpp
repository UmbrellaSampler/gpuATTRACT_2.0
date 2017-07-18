/*
 * emATTRACT.cpp
 *
 *  Created on: Jul 18, 2017
 *      Author: uwe
 */

#include "CPU_6D_EnergyService.h"
#include "emATTRACT.tpp"
#include "Chunk.h"

#ifdef CUDA
#include "GPU_6D_EnergyService.h"
#endif

namespace as {

template
class emATTRACT<CPU_6D_EnergyService<float>>;

template
class emATTRACT<CPU_6D_EnergyService<double>>;

#ifdef CUDA
template
class emATTRACT<GPU_6D_EnergyService<float>>;

template
class emATTRACT<GPU_6D_EnergyService<double>>;
#endif

}  // namespace as


