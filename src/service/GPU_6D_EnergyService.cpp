/*
 * GPU_6D_EnergyService.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include "DeviceAllocator.tpp"
#include "HostPinnedAllocator.tpp"
#include "GPU_6D_EnergyService.tpp"

namespace as {

template
class GPU_6D_EnergyService<float>;

template
class GPU_6D_EnergyService<double>;

}  // namespace as


#endif
