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
#include "DOFConverter.h"
#include "readFile.h"
#include "Server.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "SimParam.h"
#include "TypeMap.h"
#include "DataManager.h"
#include "DOFTransform.h"
#include "nativeTypesMath.h"
#include "ServiceFactory.h"

namespace as {

template<typename SERVICE>
void Configurator_6D<SERVICE>::init(CmdArgs const& args) {

	/* load dataItems */
	auto receptor = createProteinFromPDB<real_t>(args.recName);
	auto ligand = createProteinFromPDB<real_t>(args.ligName);
	auto paramTable = createParamTableFromFile<real_t>(args.paramsName);
	auto gridRec = createGridFromGridFile<real_t>(args.gridRecName);

	auto simParam = std::make_shared<SimParam<real_t>>();
	if (args.dielec == "variable") {
		simParam->dielec = Dielec::variable;
	} else {
		simParam->dielec = Dielec::constant;
	}
	simParam->epsilon = static_cast<real_t>(args.epsilon);
	simParam->ffelec = 	static_cast<real_t>(FELEC/args.epsilon);


	/* apply mapping according to receptor grid alphabet to ligand */
	auto mapVecRec = readGridAlphabetFromFile(args.alphabetRecName); // map: std::vector<unsigned>
	TypeMap typeMapRec = createTypeMapFromVector(mapVecRec);
	ligand->setNumMappedTypes(1);
	ligand->getOrCreateMappedPtr();
	applyDefaultMapping(ligand->numAtoms(), ligand->type(), ligand->type());
	applyMapping(typeMapRec, ligand->numAtoms(), ligand->type(), ligand->mappedType());


	/* read dof file */
	DOFHeader<real_t> h = readDOFHeader<real_t>(args.dofName);
	/* check file. only a receptor-ligand pair (no multi-bodies!) is allowed */
	if(!h.auto_pivot && h.pivots.size() > 2) {
		throw std::logic_error("DOF-file contains definitions for more than two molecules. Multi-body docking is not supported.");
	}


	std::vector<std::vector<DOF>> DOF_molecules_dof = readDOF(args.dofName);
	std::vector<std::vector<DOF_6D<real_t>>> DOF_molecules = DOFConverter<real_t>(DOF_molecules_dof);
	if(DOF_molecules.size() != 2) {
		throw std::logic_error("DOF-file contains definitions for more than two molecules. Multi-body docking is not supported.");
	}


	/* apply pivoting to proteins */
		if(h.auto_pivot) {
			if (!h.pivots.empty()) {
				throw std::logic_error("Auto pivot specified, but explicitly defined pivots available. (File " + args.dofName + ")" );
			}
			receptor->auto_pivotize();
			ligand->auto_pivotize();
			h.pivots.push_back(receptor->pivot());
			h.pivots.push_back(ligand->pivot());
		} else {
			if (h.pivots.size() != 2) {
				throw std::logic_error("No auto pivot specified, but number of defined pivots is incorrect. (File " + args.dofName + ")" );
			}
			receptor->pivotize(h.pivots[0]);
			ligand->pivotize(h.pivots[1]);
		}

	/* transform ligand dofs assuming that the receptor is always centered in the origin */
	transformDOF_glob2rec(DOF_molecules[0], DOF_molecules[1], h.pivots[0], h.pivots[1], h.centered_receptor, h.centered_ligands);


	/* init dof and result buffer */
	this->_dofs = std::vector<input_t>(DOF_molecules[1].size());
	for (size_t i = 0; i < DOF_molecules[1].size(); ++i) {
		this->_dofs[i].pos = DOF_molecules[1][i].pos;
		this->_dofs[i].ang = DOF_molecules[1][i].ang;
	}



	/* apply grid displacement */
	gridRec->translate(-make_real3(receptor->pivot().x,receptor->pivot().y,receptor->pivot().z));


	/* add items to dataMng */
	std::shared_ptr<DataManager> dataManager = std::make_shared<DataManager>();
	this->_ids.recId = dataManager->add(receptor);
	this->_ids.ligId = dataManager->add(ligand);
	this->_ids.gridId = dataManager->add(gridRec);
	this->_ids.tableId = dataManager->add(paramTable);
	this->_ids.paramsId = dataManager->add(simParam);

#ifdef CUDA
	if (args.deviceIds.size() > 0) {
		for (auto id : args.deviceIds) {
			dataManager->attachAllDataToDevice(id);
		}
	}
#endif

	ServiceType serviceType;
	if (args.numCPUs > 0) {
		serviceType = ServiceType::CPUEnergyService6D;
	}
#ifdef CUDA
	else {
		serviceType = ServiceType::GPUEnergyService6D;
	}
#endif

	std::shared_ptr<service_t> service = std::move(std::static_pointer_cast<service_t>(ServiceFactory::create<real_t>(serviceType, dataManager, args)));

	this->_server = std::unique_ptr<server_t>(new server_t(service));
	if (args.numCPUs > 0) {
		this->_server->createWorkers(args.numCPUs);
	} else {
		this->_server->createWorkers(args.deviceIds.size());
	}

}

template<typename SERVICE>
void Configurator_6D<SERVICE>::finalize() noexcept {

}

}  // namespace as




#endif /* SRC_CONFIGURATOR_6D_TPP_ */
