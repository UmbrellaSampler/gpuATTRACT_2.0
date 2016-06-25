/*
 * DataManager.tpp
 *
 *  Created on: May 22, 2016
 *      Author: uwe
 */

#ifndef SRC_DATAMANAGER_TPP_
#define SRC_DATAMANAGER_TPP_


#include <stdexcept>
#include "DataManager.h"
#include "DeviceManager.h"

using namespace as;

template<typename REAL>
class DataManager<REAL>::DeviceOccupancy {

	void addGrid(gridId_t id) noexcept {
		_gridIds.insert(id);
	}

	void addProtein(proteinId_t id) noexcept {
		_proteinIds.insert(id);
	}

	void addTable(tableId_t id) noexcept {
		_tableIds.insert(id);
	}

	void addParams(paramId_t id) noexcept {
		_paramIds.insert(id);
	}

	std::set<gridId_t> grids() noexcept {
		return _gridIds;
	}

	std::set<proteinId_t> proteins() noexcept {
		return _proteinIds;
	}

	std::set<tableId_t> tables() noexcept {
		return _tableIds;
	}

	std::set<paramId_t> params() noexcept {
		return _paramIds;
	}

private:
	std::set<gridId_t> _gridIds;
	std::set<proteinId_t> _proteinIds;
	std::set<tableId_t> _tableIds;
	std::set<paramId_t> _paramIds;
};

template<typename REAL>
class DataManager<REAL>::DeviceLocationMap {

	void addGrid(gridId_t id, deviceId_t dId) noexcept {
		_gridDeviceMap[id]= dId;
	}

	void addProtein(proteinId_t id, deviceId_t dId) noexcept {
		_proteinDeviceMap[id]= dId;
	}

	void addTable(tableId_t id, deviceId_t dId) noexcept {
		_tableDeviceMap[id]= dId;
	}

	void addParams(paramId_t id, deviceId_t dId) noexcept {
		_paramDeviceMap[id]= dId;
	}

	std::set<deviceId_t> idsOfGrid(gridId_t id) {
		return _gridDeviceMap.at(id);
	}

	std::set<deviceId_t> idsOfProtein(proteinId_t id) {
		return _proteinDeviceMap.at(id);
	}

	std::set<deviceId_t> idsOfTable(tableId_t id) {
		return _tableDeviceMap.at(id);
	}

	std::set<deviceId_t> idsOfParam(paramId_t id) {
		return _paramDeviceMap.at(id);
	}

private:
	std::map<gridId_t, std::set<deviceId_t>> _gridDeviceMap;
	std::map<proteinId_t, std::set<deviceId_t>> _proteinDeviceMap;
	std::map<tableId_t, std::set<deviceId_t>> _tableDeviceMap;
	std::map<paramId_t, std::set<deviceId_t>> _paramDeviceMap;
};

template<typename REAL>
class DataManager<REAL>::DataIDs {
public:

	bool operator<(DataIDs const& rhs) {
		return gridId < rhs.gridId
				|| (gridId    == rhs.gridId    && (proteinId < rhs.proteinId
				|| (proteinId == rhs.proteinId && (tableId   < rhs.tableId
				|| (tableId   == rhs.tableId   && (paramId        < rhs.paramId))))));

	}

	bool operator==(DataIDs const& rhs) {
		return gridId == rhs.gridId
				&& proteinId == rhs.proteinId
				&& tableId == rhs.tableId
				&& paramId = rhs.paramId;
	}

	gridId_t gridId;
	proteinId_t proteinId;
	tableId_t tableId;
	paramId_t paramId;
};


#define ADD_OBJECT( name ) do {                                                \
		if (name.get() == nullptr) {                                           \
			std::string obj = #name;									       \
			throw std::invalid_argument("Invalid " + obj + " (nullptr).");     \
		}                                                                      \
		_##name##s.push_back(name);                                            \
		name##Id_t id = _##name##s.size()-1;                                   \
		return id;                                                             \
	} while(0)

template<typename REAL>
proteinId_t DataManager<REAL>::addProtein(std::shared_ptr<protein_t> protein) {
	ADD_OBJECT(protein);
}

template<typename REAL>
gridId_t DataManager<REAL>::addGridUnion(std::shared_ptr<gridUnion_t> grid) {
	ADD_OBJECT(grid);
}

template<typename REAL>
tableId_t DataManager<REAL>::addParamTable(std::shared_ptr<table_t> table) {
	ADD_OBJECT(table);
}

template<typename REAL>
paramId_t DataManager<REAL>::addSimParam(std::shared_ptr<param_t> param) {
	ADD_OBJECT(param);
}

#define CHECK_VALID_ID( name ) do {                                           \
		if(name##Id < _##name##s.size() && _##name##s[name##Id] != nullptr) {  \
			std::string obj = #name;                                      \
			std::stringstream stream; \
			stream << name##Id; \
			throw std::invalid_argument("Invalid " + obj + "Id (" + stream.str() + ").");        \
		}                                                                 \
	} while(0)

#define GET_OBJECT( name )  do { \
		CHECK_VALID_ID(name); \
		return _##name##s[name##Id]; \
	} while(0)

template<typename REAL>
auto DataManager<REAL>::getProtein(proteinId_t proteinId) -> std::shared_ptr<protein_t>{
	GET_OBJECT(protein);
}

template<typename REAL>
auto DataManager<REAL>::getGridUnion(gridId_t gridId) -> std::shared_ptr<gridUnion_t> {
	GET_OBJECT(grid);
}

template<typename REAL>
auto DataManager<REAL>::getParamTable(tableId_t tableId) -> std::shared_ptr<table_t> {
	GET_OBJECT(table);
}

template<typename REAL>
auto DataManager<REAL>::getSimParam(paramId_t paramId) -> std::shared_ptr<param_t> {
	GET_OBJECT(param);
}

#define REMOVE_OBJECT(name) do { \
	CHECK_VALID_ID(name); \
	_##name##s[name##Id] = nullptr; \
} while(0)

template<typename REAL>
void DataManager<REAL>::removeProtein(proteinId_t proteinId) {
	REMOVE_OBJECT(protein);
}

template<typename REAL>
void DataManager<REAL>::removeGridUnion(gridId_t gridId) {
	REMOVE_OBJECT(grid);
}

template<typename REAL>
void DataManager<REAL>::removeParamTable(tableId_t tableId) {
	REMOVE_OBJECT(table);
}

template<typename REAL>
void DataManager<REAL>::removeSimParam(paramId_t paramId) {
	REMOVE_OBJECT(param);
}

template<typename REAL>
void DataManager<REAL>::removeEnsemble(std::shared_ptr<ensemble_t>& ensemble) {
	auto proteinIds = ensemble->ids();
	for (auto id : proteinIds) {
		removeProtein(id);
	}
}

template<typename REAL>
void DataManager<REAL>::removeAll() {
	_grids.clear();
	_proteins.clear();
	_tables.clear();
	_params.clear();
}

#ifdef CUDA
#include <cuda_runtime.h>

template<typename REAL>
void DataManager<REAL>::attachGridUnionToDevice(gridId_t gridId, deviceId_t deviceId) {
	CHECK_VALID_ID(grid);
	checkDeviceIdAndSet(deviceId);
	//TODO: finish implementation
}

template<typename REAL>
void DataManager<REAL>::detachGridUnionFromDevice(gridId_t gridId, deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::attachProteinToDevice(proteinId_t proteinId, deviceId_t deviceId) {
	CHECK_VALID_ID(protein);
	checkDeviceIdAndSet(deviceId);
	auto& sh_protein = _proteins[proteinId];
	DeviceManager<REAL>::createDeviceProtein(sh_protein.get());
}

template<typename REAL>
void DataManager<REAL>::detachProteinFromDevice(proteinId_t proteinId, deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::attachParamTableToDevice(tableId_t tableId, deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::detachParamTableFromDevice(tableId_t paramId, deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::attachSimParamToDevice(paramId_t paramId, deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::detachSimParamFromDevice(paramId_t paramId, deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::attachAllDataToDevice(deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::attachDataToAllDevices() {}

template<typename REAL>
void DataManager<REAL>::releaseDevice(deviceId_t deviceId) {}

template<typename REAL>
void DataManager<REAL>::releaseAllDevices() {}

template<typename REAL>
void DataManager<REAL>::checkDeviceIdAndSet(deviceId_t deviceId) {
	cudaError_t err = cudaSetDevice(deviceId);
	if(err!=cudaSuccess) {
		std::stringstream stream;
		stream << deviceId;
		throw std::invalid_argument("Invalid deviceId (" + stream.str() + ").");
	}
}
#endif





#endif /* SRC_DATAMANAGER_TPP_ */
