/*
 * GPU_6D_EnergyService.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include "DeviceAllocator.tpp"
#include "HostPinnedAllocator.tpp"
#include "GPUEnergyService6D.tpp"

namespace as {

template
class GPUEnergyService6D<float>;

template
class GPUEnergyService6D<double>;

}  // namespace as


#endif
