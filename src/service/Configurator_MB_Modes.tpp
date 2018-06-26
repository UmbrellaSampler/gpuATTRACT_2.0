/*
 * Configurator_6D_Modes.tpp
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATOR_MB_MODES_TPP_
#define SRC_SERVICE_CONFIGURATOR_MB_MODES_TPP_

#include <exception>
#include <vector>

#include "Configurator_MB_Modes.h"
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
#include "scoring_kernel.h"
#include <iostream>
namespace as {

template<typename SERVICE>
void Configurator_MB_Modes<SERVICE>::init(CmdArgs const& args) {
	std::shared_ptr<DataManager> dataManager = std::make_shared<DataManager>();


	/* read dof file */
	DOFHeader<real_t> h = readDOFHeader<real_t>(args.dofName);

	/* check file. only a receptor-ligand pair (no multi-bodies!) is allowed */
	if(!h.auto_pivot && h.pivots.size() > 2) {
		throw std::logic_error("DOF-file contains definitions for more than two molecules. Multi-body docking is not supported.");
	}
	// TODO: transform DOF_6D to input_t
	//std::vector<std::vector<DOF_6D_Modes<real_t>>> DOF_molecules = std::vector<std::vector<DOF_6D_Modes<real_t>>>();
	std::vector<std::vector<DOF>> DOF_molecules_dof = readDOF(args.dofName);
	//std::vector<DOF_MB_Modes<real_t>> DOF_molecules = DOFConverter_MB_Modes<real_t>(DOF_molecules_dof);
	this->_dofs = DOFConverter_MB_Modes<real_t>(DOF_molecules_dof);

	/* load dataItems */
	this->_ids.radius_cutoff = args.radius_cutoff;
	this->_ids.numProteins = args.numProtein;
	auto mapVec = readGridAlphabetFromFile(args.alphabetName); // map: std::vector<unsigned>
	TypeMap typeMap = createTypeMapFromVector(mapVec);
	int numProtein = 2;
	for( int idxProtein = 0; idxProtein < numProtein; ++idxProtein){
		std::cout <<args.proteinNames[idxProtein] << std::endl;
		auto protein = createProteinFromPDB<real_t>(args.proteinNames[idxProtein]);
		auto grid = createGridFromGridFile<real_t>(args.gridNames[idxProtein]);
		protein->setNumModes( args.numModes );
		protein->auto_pivotize();
		protein->setNumMappedTypes(1);
		readHMMode<real_t>(protein, args.modeNames[idxProtein]);
		protein->getOrCreateMappedPtr();
		protein-> scaleModeEigenValues( args.modeEVFactor );
		applyDefaultMapping(protein->numAtoms(), protein->type(), protein->type());
		applyMapping(typeMap, protein->numAtoms(), protein->type(), protein->mappedType());
		Common_MB_Modes::numModes[idxProtein] = args.numModes;
		auto id_grid = dataManager->add(grid);
		auto id_protein = dataManager->add(protein);
		ProtConfig config( id_grid , id_protein, idxProtein, false, protein->pivot());
		this->_ids.proteins.push_back( config );
	}

		if(h.auto_pivot) {
			if (!h.pivots.empty()) {
				throw std::logic_error("Auto pivot specified, but explicitly defined pivots available. (File " + args.dofName + ")" );
			}
			//receptor->auto_pivotize();
			//ligand->auto_pivotize();
			//h.pivots.push_back(receptor->pivot());
			//h.pivots.push_back(ligand->pivot());
		} else {
			if (h.pivots.size() != 2) {
				throw std::logic_error("No auto pivot specified, but number of defined pivots is incorrect. (File " + args.dofName + ")" );
			}
			//receptor->pivotize(h.pivots[0]);
			//ligand->pivotize(h.pivots[1]);
		}


		auto paramTable = createParamTableFromFile<real_t>(args.paramsName);

		auto simParam = std::make_shared<SimParam<real_t>>();
			if (args.dielec == "variable") {
				simParam->dielec = Dielec::variable;
			} else {
				simParam->dielec = Dielec::constant;
			}
			simParam->epsilon = static_cast<real_t>(args.epsilon);
			simParam->ffelec = 	static_cast<real_t>(FELEC/args.epsilon);


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
		serviceType = ServiceType::CPUEnergyServiceMBModes;
	}

	// TODO: ServiceType::GPUEnergyService6DModes is not yet available
#ifdef CUDA
	else {
		serviceType = ServiceType::GPUEnergyServiceMBModes;
	}
#endif

	std::shared_ptr<service_t> service = std::move(std::static_pointer_cast<service_t>(ServiceFactory::create<real_t>(serviceType, dataManager, args)));

	this->_server = std::unique_ptr<server_t>(new server_t(service));
	if (args.numCPUs > 0) {
		this->_server->createWorkers(args.numCPUs);
	} else {
		//this->_server->createWorkers(args.deviceIds.size());
		this->_server->createWorkers(1);
	}

}

template<typename SERVICE>
void Configurator_MB_Modes<SERVICE>::finalize() noexcept {

}

}  // namespace as



#endif /* SRC_SERVICE_CONFIGURATOR_6D_MODES_TPP_ */
