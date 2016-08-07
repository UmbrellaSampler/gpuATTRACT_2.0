/*
 * DeviceDataConfigurator.cpp
 *
 *  Created on: Jul 19, 2016
 *      Author: uwe
 */

#ifdef CUDA

#include <stdexcept>
#include <sstream>
#include "cuda_runtime.h"
#include "DeviceDataConfigurator.h"

namespace as {

void DeviceDataConfigurator::checkDeviceIdAndSetCurrent(deviceId_t deviceId) {
	cudaError_t err = cudaSetDevice(deviceId);
	if(err!=cudaSuccess) {
		std::stringstream stream;
		stream << deviceId;
		std::string(cudaGetErrorString(err));
		throw std::invalid_argument("CUDA error: Invalid deviceId (" + stream.str() + ").");
	}
}

}

#include "DeviceDataConfigurator.tpp"

namespace as {

template
std::shared_ptr<DeviceProtein<float>> DeviceDataConfigurator::attach(const std::shared_ptr<Protein<float>> protein, deviceId_t id);

template
std::shared_ptr<DeviceProtein<double>> DeviceDataConfigurator::attach(const std::shared_ptr<Protein<double>> protein, deviceId_t id);

template
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceProtein<float>>, deviceId_t);

template
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceProtein<double>>, deviceId_t);

template
std::shared_ptr<DeviceGridUnion<float>> DeviceDataConfigurator::attach(const std::shared_ptr<GridUnion<float>> protein, deviceId_t id);

template
std::shared_ptr<DeviceGridUnion<double>> DeviceDataConfigurator::attach(const std::shared_ptr<GridUnion<double>> protein, deviceId_t id);

template
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceGridUnion<float>>, deviceId_t);

template
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceGridUnion<double>>, deviceId_t);

}

#endif
