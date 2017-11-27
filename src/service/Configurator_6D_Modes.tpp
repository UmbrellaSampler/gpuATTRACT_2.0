#ifndef SRC_CONFIGURATOR_6D_MODES_TPP_
#define SRC_CONFIGURATOR_6D_MODES_TPP_

#include <exception>
#include <vector>

#include "Configurator_6D_Modes.h"

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
#include "ServiceFactory.h"


namespace as {

template<typename SERVICE>
void Configurator_6D_Modes<SERVICE>::init(CmdArgs const& args) noexcept {

	/* load dataItems */
	auto receptor = createProteinFromPDB<real_t>(args.recName);
	auto ligand = createProteinFromPDB<real_t>(args.ligName);
	auto gridRec = createGridFromGridFile<real_t>(args.gridNameRec);
	auto gridLig = createGridFromGridFile<real_t>(args.gridNameLig);

	auto paramTable = createParamTableFromFile<real_t>(args.paramsName);

	if(args.numModes > 0){
		receptor->setNumModes(args.numModes);
		ligand->setNumModes(args.numModes);
		readHMMode<real_t>(receptor, args.recModesName);
		readHMMode<real_t>(ligand, args.ligModesName);
	}


	auto simParam = std::make_shared<SimParam<real_t>>();
	if (args.dielec == "variable") {
		simParam->dielec = Dielec::variable;
	} else {
		simParam->dielec = Dielec::constant;
	}
	simParam->epsilon = static_cast<real_t>(args.epsilon);
	simParam->ffelec = 	static_cast<real_t>(FELEC/args.epsilon);


	/* apply mapping according to receptor grid alphabet to ligand */
	auto mapVec = readGridAlphabetFromFile(args.alphabetNameRec); // map: std::vector<unsigned>
	TypeMap typeMap = createTypeMapFromVector(mapVec);
	ligand->setNumMappedTypes(1);
	ligand->getOrCreateMappedPtr();
	applyDefaultMapping(ligand->numAtoms(), ligand->type(), ligand->type());
	applyMapping(typeMap, ligand->numAtoms(), ligand->type(), ligand->mappedType());


	auto mapVecRec = readGridAlphabetFromFile(args.alphabetNameLig); // map: std::vector<unsigned>
	TypeMap typeMapRec = createTypeMapFromVector(mapVecRec);
	receptor->setNumMappedTypes(1);
	receptor->getOrCreateMappedPtr();
	applyDefaultMapping(receptor->numAtoms(), receptor->type(), receptor->type());
	applyMapping(typeMapRec, receptor->numAtoms(), receptor->type(), receptor->mappedType());

	/* read dof file */
	DOFHeader<real_t> h = readDOFHeader<real_t>(args.dofName);
	/* check file. only a receptor-ligand pair (no multi-bodies!) is allowed */
	if(!h.auto_pivot && h.pivots.size() > 2) {
		throw std::logic_error("DOF-file contains definitions for more than two molecules. Multi-body docking is not supported.");
	}

	std::vector<std::vector<DOF_6D_Modes<real_t>>> DOF_molecules = readDOF_6D_Modes<real_t>(args.dofName,args.numModes);
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
	_dofs=std::vector<input_t>(DOF_molecules[1].size());


	for (size_t i = 0; i < DOF_molecules[1].size(); ++i) {
		_dofs[i].dof.push_back(DOF_molecules[0][i]);
		_dofs[i].dof.push_back(DOF_molecules[1][i]);
	}

	/* apply grid displacement */
	gridRec->translate(-make_real3(h.pivots[0].x,h.pivots[0].y,h.pivots[0].z));
	gridLig->translate(-make_real3(h.pivots[1].x,h.pivots[1].y,h.pivots[1].z));

	/* add items to dataMng */
	std::shared_ptr<DataManager> dataManager = std::make_shared<DataManager>();
	_ids.recId = dataManager->add(receptor);
	_ids.ligId = dataManager->add(ligand);
	_ids.gridRecId = dataManager->add(gridRec);
	_ids.gridLigId = dataManager->add(gridLig);
	_ids.tableId = dataManager->add(paramTable);
	_ids.paramsId = dataManager->add(simParam);





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
#ifdef CUDA
	else {
		serviceType = ServiceType::GPUEnergyService6DModes;
	}
#endif




	/* init server and service*/
//	std::shared_ptr<Service<dof_t, common_t, result_t>> service;
//	if (std::is_same<service_t, CPU_6D_EnergyService<real_t>>::value) {
//		service = std::make_shared<CPU_6D_EnergyService<real_t>>();
//	}
//#ifdef CUDA
//	else if (std::is_same<service_t, GPU_6D_EnergyService<real_t>>::value){
//		 service = std::make_shared<GPU_6D_EnergyService<real_t>>(args.deviceIds);
//	}
//#endif


	//ToDo: create ServiceFactory with constructor for CmdArgs
	std::shared_ptr<service_t> service = std::move(std::static_pointer_cast<service_t>(ServiceFactory::create<real_t>(serviceType, dataManager, args)));

	service->initAllocators();
	_server = std::unique_ptr<server_t>(new server_t(service));
	if (args.numCPUs > 0) {
		_server->createWorkers(args.numCPUs);
	} else {
		_server->createWorkers(args.deviceIds.size());
	}

}

template<typename SERVICE>
void Configurator_6D_Modes<SERVICE>::finalize() noexcept {

}

}  // namespace as




#endif /* SRC_CONFIGURATOR_6D_TPP_ */
