/*
 * DataManager.h
 *
 *  Created on: Mar 30, 2016
 *      Author: uwe
 */

#ifndef SRC_DATAMANAGER_H_
#define SRC_DATAMANAGER_H_

#include <set>
#include <map>
#include <unordered_map>

//TODO: remove
#define CUDA

namespace as {

constexpr unsigned maxNumProteins = 1000;
constexpr unsigned maxNumGrids = 50;

using proteinId_t = unsigned;
using gridId_t = unsigned;
using tableId_t = unsigned;
using paramId_t = unsigned;
using deviceId_t = unsigned;

template<typename REAL>
class GridUnion;

template<typename REAL>
class Protein;

template<typename REAL>
class ParamTable;

template<typename REAL>
class SimParam;

template<typename REAL>
class Ensemble;

template<typename REAL>
class DataManager {

private:
	using gridUnion_t = GridUnion<REAL>;
	using protein_t = Protein<REAL>;
	using table_t = ParamTable<REAL>;
	using param_t = SimParam<REAL>;
	using ensemble_t = Ensemble<REAL>;

	class DeviceOccupancy;
	class DeviceLocationMap;
	class DataIDs;

public:
	DataManager() {};
	~DataManager() {};

	proteinId_t addProtein(std::shared_ptr<protein_t> protein);
	gridId_t addGridUnion(std::shared_ptr<gridUnion_t> gridUnion);
	tableId_t addParamTable(std::shared_ptr<table_t> table);
	paramId_t addSimParam(std::shared_ptr<param_t> simPar);
	std::vector<proteinId_t> addEnsemble(std::shared_ptr<ensemble_t>& ensemble);

	std::shared_ptr<protein_t> getProtein(proteinId_t proteinId);
	std::shared_ptr<gridUnion_t> getGridUnion(gridId_t gridId);
	std::shared_ptr<table_t> getParamTable(tableId_t tableId);
	std::shared_ptr<param_t> getSimParam(paramId_t paramId);

	void removeProtein(proteinId_t proteinId);
	void removeGridUnion(gridId_t gridId);
	void removeParamTable(tableId_t tableId);
	void removeSimParam(paramId_t paramId);
	void removeEnsemble(std::shared_ptr<ensemble_t>& ensemble);

	void removeAll();

#ifdef CUDA
	void attachGridUnionToDevice(gridId_t gridId, deviceId_t deviceId);
	void detachGridUnionFromDevice(gridId_t gridId, deviceId_t deviceId);

	void attachProteinToDevice(proteinId_t proteinId, deviceId_t deviceId);
	void detachProteinFromDevice(proteinId_t proteinId, deviceId_t deviceId);

	/* Parameter Objects: parameter table (AttrParamTable)
	 * and simulation parameter (simParam) */
	void attachParamTableToDevice(tableId_t tableId, deviceId_t deviceId);
	void detachParamTableFromDevice(tableId_t paramId, deviceId_t deviceId);

	void attachSimParamToDevice(paramId_t paramId, deviceId_t deviceId);
	void detachSimParamFromDevice(paramId_t paramId, deviceId_t deviceId);

	void attachAllDataToDevice(deviceId_t deviceId);
	void attachDataToAllDevices();

	void releaseDevice(deviceId_t deviceId);
	void releaseAllDevices();
#endif

private:

	std::vector<std::shared_ptr<protein_t>> _proteins;
	std::vector<std::shared_ptr<gridUnion_t>> _grids;

	std::vector<std::shared_ptr<table_t>> _tables;
	std::vector<std::shared_ptr<param_t>> _params;

#ifdef CUDA
	/** Maps deviceId to a set of objects (grids, proteins, ...) */
	std::map<deviceId_t, DeviceOccupancy> _deviceOccupancyMap;

	/** Maps a objectId (gridId_t, proteinId_t, ...) to a set of devices the object is attached to */
	DeviceLocationMap _deviceLocationMap;


	/** Lookup table for each dataID-collection (grid, protein, table, param).
	 * The look-up returns a set of device IDs
	 * indicating the devices that can be used to process collection.*/
	std::map<DataIDs, std::set<unsigned> > _commonDeviceID_LookUp;

	/**
	 * Checks if a cuda device with deviceId is available. If yes, this device is
	 * set to the active device. Otherwise, an exception is thrown.
	 */
	void checkDeviceIdAndSet(deviceId_t deviceId);

#endif

};

}  // namespace as



#endif /* SRC_DATAMANAGER_H_ */
