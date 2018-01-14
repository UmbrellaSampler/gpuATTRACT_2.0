/*
 * ConfiguratorBase.tpp
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATORBASE_TPP_
#define SRC_SERVICE_CONFIGURATORBASE_TPP_

#include "ConfiguratorBase.h"
#include "SimulationInstances.h"
#include "readFile.h"

namespace as {

template<typename GenericTypes>
auto ConfiguratorBase<GenericTypes>::getBaseSimulationAttributes(CmdArgs const& args) -> SimulationConfiguration<real_t> {
	SimulationConfiguration<real_t> instances;
	instances.proteins.add(createProteinFromPDB<real_t>(args.recName));
	instances.proteins.add(createProteinFromPDB<real_t>(args.ligName));
	instances.paramTable = createParamTableFromFile<real_t>(args.paramsName);
	instances.grids.add(createGridFromGridFile<real_t>(args.gridRecName));
	instances.proteinAlphabets.add(readGridAlphabetFromFile(args.alphabetRecName));
	instances.simParam = std::make_shared<SimParam<real_t>>();

	if (args.dielec == "variable") {
		instances.simParam->dielec = Dielec::variable;
	} else {
		instances.simParam->dielec = Dielec::constant;
	}
	instances.simParam->epsilon = static_cast<real_t>(args.epsilon);
	instances.simParam->ffelec = static_cast<real_t>(FELEC/args.epsilon);


}

} // namespace as



#endif /* SRC_SERVICE_CONFIGURATORBASE_TPP_ */
