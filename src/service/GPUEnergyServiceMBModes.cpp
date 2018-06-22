/*
 * GPU_6D_EnergyService.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include "DeviceAllocator.tpp"
#include "HostPinnedAllocator.tpp"
#include "GPUEnergyServiceMBModes.tpp"

namespace as {

template
class GPUEnergyServiceMBModes<float>;

template
class GPUEnergyServiceMBModes<double>;

}  // namespace as


#endif
