/*
 * CPUEnergyService6D_MB_Modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef CPUENERGYSERVICE6D_MB_MODES_H_
#define CPUENERGYSERVICE6D_MB_MODES_H_


#include "CPUEnergyService.h"
#include "publicTypes.h"
#include "Types_6D_MB_Modes.h"
#include "nativeTypesWrapper.h"

namespace as {

template<typename SERVICE>
class Configurator_6D_MB_Modes;

template<typename REAL>
class CPUEnergyService6D_MB_Modes : public CPUEnergyService<Types_6D_MB_Modes<REAL>> {
public:
	using typename CPUEnergyService<Types_6D_MB_Modes<REAL>>::service_t;
	using real_t = typename TypeWrapper<REAL>::real_t;
	using typename service_t::input_t;
	using typename service_t::common_t;
	using typename service_t::result_t;
	using configurator_t = Configurator_6D_MB_Modes<real_t>;

	/* need public for instantiate the server in configurator */

	using typename service_t::workItem_t;
	using typename service_t::itemProcessor_t;

	explicit CPUEnergyService6D_MB_Modes(std::shared_ptr<DataManager> dataMng);
	virtual ~CPUEnergyService6D_MB_Modes() {};

	itemProcessor_t createItemProcessor() override;

private:

	class Buffer;

};

}  // namespace as





#endif /* CPUENERGYSERVICE6DMODES_H_ */
