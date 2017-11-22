/*
 * ServiceType.h
 *
 *  Created on: Aug 12, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_SERVICETYPE_H_
#define SRC_SERVICE_SERVICETYPE_H_

namespace as {

enum class ServiceType {
	CPUEnergyService6D,
	CPUEnergyService6DModes,
#ifdef CUDA
	GPUEnergyService6D,
	GPUEnergyService6DModes
#endif
};

} // namespace



#endif /* SRC_SERVICE_SERVICETYPE_H_ */
