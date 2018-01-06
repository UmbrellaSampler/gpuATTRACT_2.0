/*
 * GPU_6D_MB_EnergyService.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include "DeviceAllocator.tpp"
#include "HostPinnedAllocator.tpp"
#include "GPUEnergyService6D_MB_Modes.tpp"

namespace as {

template
class GPUEnergyService6D_MB_Modes<float>;

template
class GPUEnergyService6D_MB_Modes<double>;

}  // namespace as


#endif
