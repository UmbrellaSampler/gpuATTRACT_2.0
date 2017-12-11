/*
 * GPU_6D_EnergyService.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include "DeviceAllocator.tpp"
#include "HostPinnedAllocator.tpp"
#include "GPUEnergyService6DModes.tpp"

namespace as {

template
class GPUEnergyService6DModes<float>;

template
class GPUEnergyService6DModes<double>;

}  // namespace as


#endif
