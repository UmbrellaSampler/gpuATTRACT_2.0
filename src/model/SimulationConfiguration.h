/*
 * SimulationAttributes.h
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVER_SIMULATIONCONFIGURATION_H_
#define SRC_SERVER_SIMULATIONCONFIGURATION_H_

#include <vector>
#include <memory>
#include "Protein.h"
#include "ParamTable.h"
#include "SimParam.h"
#include "DOFHeader.h"
#include "TypeMap.h"

namespace as {

template<typename REAL>
struct SimulationConfiguration {
	std::vector<std::shared_ptr<Protein<REAL>>> proteins;
	std::vector<std::shared_ptr<Protein<REAL>>> grids;
	std::vector<TypeMap> proteinAlphabets;
	std::vector<std::shared_ptr<ParamTable<REAL>>> paramTable;
	std::shared_ptr<SimParam<REAL>> simParam;
	DOFHeader<REAL> dofHeader;
};

} // namespace as



#endif /* SRC_SERVER_SIMULATIONCONFIGURATION_H_ */
