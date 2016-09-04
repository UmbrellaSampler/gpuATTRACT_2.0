/*
 * Configurator_6D.tpp
 *
 *  Created on: Aug 17, 2016
 *      Author: uwe
 */

#ifndef SRC_CONFIGURATOR_6D_TPP_
#define SRC_CONFIGURATOR_6D_TPP_

#include <exception>
#include <vector>

#include "Configurator_6D.h"
#include "readFile.h"
#include "DOF_6D.h"
//#include "ServerIncludes.h"
#include "Server.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "SimParam.h"
#include "TypeMap.h"
#include "DataManager.h"
#include "DOFTransform.h"
#include "nativeTypesMath.h"

namespace as {

template<typename REAL>
void Configurator_6D<REAL>::init(CmdArgs const& args) {

	/* load dataItems */
	auto receptor = createProteinFromPDB<REAL>(args.recName);
	auto ligand = createProteinFromPDB<REAL>(args.ligName);
	auto grid = createGridFromGridFile<REAL>(args.gridName);
	auto paramTable = createParamTableFromFile<REAL>(args.paramsName);

	auto simParam = std::make_shared<SimParam<REAL>>();
	if (args.dielec == "variable") {
		simParam->dielec = Dielec::variable;
	} else {
		simParam->dielec = Dielec::constant;
	}
	simParam->epsilon = static_cast<REAL>(args.epsilon);
	simParam->ffelec = 	static_cast<REAL>(FELEC/args.epsilon);


	/* apply mapping according to receptor grid alphabet to ligand */
	auto mapVec = readGridAlphabetFromFile(args.alphabetName); // map: std::vector<unsigned>
	TypeMap typeMap = createTypeMapFromVector(mapVec);
	ligand->setNumMappedTypes(1);
	ligand->getOrCreateMappedPtr();
	applyDefaultMapping(ligand->numAtoms(), ligand->type(), ligand->type());
	applyMapping(typeMap, ligand->numAtoms(), ligand->type(), ligand->mappedType());

	/* add items to dataMng */
	std::shared_ptr<DataManager> dataManager = std::make_shared<DataManager>();
	_ids.recId = dataManager->add(receptor);
	_ids.ligId = dataManager->add(ligand);
	_ids.gridId = dataManager->add(grid);
	_ids.tableId = dataManager->add(paramTable);
	_ids.paramsId = dataManager->add(simParam);

#ifdef CUDA
	if (args.deviceIds.size() > 0) {
		for (auto id : args.deviceIds) {
			dataManager->attachAllDataToDevice(id);
		}
	}
#endif

	/* read dof file */
	DOFHeader<REAL> h = readDOFHeader<REAL>(args.dofName);
	/* check file. only a receptor-ligand pair (no multi-bodies!) is allowed */
	if(!h.auto_pivot && h.pivots.size() > 2) {
		throw std::logic_error("DOF-file contains definitions for more than two molecules. Multi-body docking is not supported.");
	}

	std::vector<std::vector<DOF_6D<REAL>>> DOF_molecules = readDOF_6D<REAL>(args.dofName);
	if(DOF_molecules.size() != 2) {
		throw std::logic_error("DOF-file contains definitions for more than two molecules. Multi-body docking is not supported.");
	}

	/* transform ligand dofs assuming that the receptor is always centered in the origin */
	transformDOF_glob2rec(DOF_molecules[0], DOF_molecules[1], h.pivots[0], h.pivots[1], h.centered_receptor, h.centered_ligands);


	/* init dof and result buffer */
	_dofs = std::vector<dof_t>(DOF_molecules[1].size());
	for (size_t i = 0; i < DOF_molecules[1].size(); ++i) {
		_dofs[i].pos = DOF_molecules[1][i].pos;
		_dofs[i].ang = DOF_molecules[1][i].ang;
	}

	/* init server and service*/
	std::shared_ptr<service_t> service = std::make_shared<service_t>();
	service->initAllocators();
	service->setDataManager(dataManager); // do not forget!
	_server = std::unique_ptr<server_t>(new server_t(service));
	if (args.numCPUs > 0) {
		_server->createWorkers(args.numCPUs);
	} else {
		_server->createWorkers(args.deviceIds.size());
	}

	/* apply pivoting to proteins */
	if(h.auto_pivot) {
		if (!h.pivots.empty()) {
			throw std::logic_error("Auto pivot specified, but explicitly defined pivots available. (File " + args.dofName + ")" );
		}
		h.pivots.push_back(receptor->pivot());
		h.pivots.push_back(ligand->pivot());
	} else {
		if (h.pivots.size() != 2) {
			throw std::logic_error("No auto pivot specified, but number of defined pivots is incorrect. (File " + args.dofName + ")" );
		}
		receptor->pivotize(h.pivots[0]);
		ligand->pivotize(h.pivots[1]);
	}

	/* apply grid displacement */
	auto pivot = h.pivots[0];
	grid->translate(-make_real3(pivot.x,pivot.y,pivot.z));

}

template<typename REAL>
void Configurator_6D<REAL>::finalize() {

}

}  // namespace as




#endif /* SRC_CONFIGURATOR_6D_TPP_ */
