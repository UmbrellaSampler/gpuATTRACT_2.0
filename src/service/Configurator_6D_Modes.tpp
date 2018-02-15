/*
 * Configurator_6D_Modes.tpp
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATOR_6D_MODES_TPP_
#define SRC_SERVICE_CONFIGURATOR_6D_MODES_TPP_

#include <exception>
#include <vector>

#include "Configurator_6D_Modes.h"
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
void Configurator_6D_Modes<SERVICE>::init(CmdArgs const& args) {

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

	// TODO: transform DOF_6D to input_t
	//std::vector<std::vector<DOF_6D_Modes<real_t>>> DOF_molecules = std::vector<std::vector<DOF_6D_Modes<real_t>>>();
	std::vector<std::vector<DOF>> DOF_molecules_dof = readDOF(args.dofName);
	std::vector<std::vector<DOF_6D_Modes<real_t>>> DOF_molecules = DOFConverter_Modes<real_t>(DOF_molecules_dof);
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
		this->_dofs[i]._6D.pos = DOF_molecules[1][i]._6D.pos;
		this->_dofs[i]._6D.ang = DOF_molecules[1][i]._6D.ang;
		std::copy(DOF_molecules[1][i].modesRec, DOF_molecules[1][i].modesRec + receptor->numModes(), this->_dofs[i].modesRec);
		std::copy(DOF_molecules[1][i].modesLig, DOF_molecules[1][i].modesLig + ligand->numModes(), this->_dofs[i].modesLig);
	}



	/* apply grid displacement */
	gridRec->translate(-make_real3(receptor->pivot().x,receptor->pivot().y,receptor->pivot().z));


	/* add items to dataMng */
	std::shared_ptr<DataManager> dataManager = std::make_shared<DataManager>();
	this->_ids.recId = dataManager->add(receptor);
	this->_ids.ligId = dataManager->add(ligand);
	this->_ids.gridIdRec = dataManager->add(gridRec);
	this->_ids.tableId = dataManager->add(paramTable);
	this->_ids.paramsId = dataManager->add(simParam);



	receptor->setNumModes(args.numModes);
	ligand->setNumModes(args.numModes);


	readHMMode<real_t>(receptor, args.recModesName);
	readHMMode<real_t>(ligand, args.ligModesName);

	auto mapVecLig = readGridAlphabetFromFile(args.alphabetLigName); // map: std::vector<unsigned>
	TypeMap typeMapLig = createTypeMapFromVector(mapVecLig);
	receptor->setNumMappedTypes(1);
	receptor->getOrCreateMappedPtr();
	applyDefaultMapping(receptor->numAtoms(), receptor->type(), receptor->type());
	applyMapping(typeMapLig, receptor->numAtoms(), receptor->type(), receptor->mappedType());

	auto gridLig = createGridFromGridFile<real_t>(args.gridLigName);
	gridLig->translate(-make_real3(ligand->pivot().x,ligand->pivot().y,ligand->pivot().z));
	this->_ids.gridIdLig = dataManager->add(gridLig);


#ifdef CUDA
	if (args.deviceIds.size() > 0) {
		for (auto id : args.deviceIds) {
			dataManager->attachAllDataToDevice(id);
		}
	}
#endif

	ServiceType serviceType;
	if (args.numCPUs > 0) {
		serviceType = ServiceType::CPUEnergyService6DModes;
	}

	// TODO: ServiceType::GPUEnergyService6DModes is not yet available
#ifdef CUDA
	else {
		serviceType = ServiceType::GPUEnergyService6DModes;
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
void Configurator_6D_Modes<SERVICE>::finalize() noexcept {

}

}  // namespace as



#endif /* SRC_SERVICE_CONFIGURATOR_6D_MODES_TPP_ */
