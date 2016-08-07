/*
 * DeviceDataConfigurator.tpp
 *
 *  Created on: Jul 19, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEDATACONFIGURATOR_TPP_
#define SRC_DEVICEDATACONFIGURATOR_TPP_

#ifdef CUDA
#include <cstring>

#include "DeviceDataConfigurator.h"
#include "cuda_runtime.h"

#include "Protein.h"
#include "DeviceProtein.h"
#include "GridUnion.h"
#include "DeviceGridUnion.h"
#include "ParamTable.h"
#include "DeviceParamTable.h"
#include "macros.h"

namespace as {

template <typename REAL>
std::shared_ptr<DeviceProtein<REAL>> DeviceDataConfigurator::attach(const std::shared_ptr<Protein<REAL>> protein, deviceId_t deviceId) {
	ASSERT(deviceId >= 0);
	checkDeviceIdAndSetCurrent(deviceId);

	unsigned numAtoms = protein->numAtoms();

	REAL *d_xPos;
	CUDA_CHECK(cudaMalloc((void**) &d_xPos, numAtoms * sizeof(REAL)));
	CUDA_CHECK(cudaMemcpy(d_xPos, protein->xPos(), numAtoms * sizeof(REAL), cudaMemcpyHostToDevice));
	REAL *d_yPos;
	CUDA_CHECK(cudaMalloc((void**) &d_yPos, numAtoms * sizeof(REAL)));
	CUDA_CHECK(cudaMemcpy(d_yPos, protein->yPos(), numAtoms * sizeof(REAL), cudaMemcpyHostToDevice));
	REAL *d_zPos;
	CUDA_CHECK(cudaMalloc((void**) &d_zPos, numAtoms * sizeof(REAL)));
	CUDA_CHECK(cudaMemcpy(d_zPos, protein->zPos(), numAtoms * sizeof(REAL), cudaMemcpyHostToDevice));

	unsigned* d_type;
	CUDA_CHECK(cudaMalloc((void**) &d_type, numAtoms * sizeof(unsigned)));
	CUDA_CHECK(cudaMemcpy(d_type, protein->type(), numAtoms * sizeof(unsigned), cudaMemcpyHostToDevice));
	unsigned* d_mappedType;
	CUDA_CHECK(cudaMalloc((void**) &d_mappedType, numAtoms * sizeof(unsigned)));
	CUDA_CHECK(cudaMemcpy(d_mappedType, protein->mappedType(), numAtoms * sizeof(unsigned), cudaMemcpyHostToDevice));
	REAL* d_charge;
	CUDA_CHECK(cudaMalloc((void**) &d_charge, numAtoms * sizeof(REAL)));
	CUDA_CHECK(cudaMemcpy(d_charge, protein->charge(), numAtoms * sizeof(REAL), cudaMemcpyHostToDevice));

	unsigned numModes = protein->numModes();
	REAL* d_xModes = NULL;
	REAL* d_yModes = NULL;
	REAL* d_zModes = NULL;
	if (numModes != 0) {
		CUDA_CHECK(cudaMalloc((void**) d_xModes, numAtoms * numModes * sizeof(REAL)));
		CUDA_CHECK(cudaMemcpy(d_xModes, protein->xModes(), numAtoms * numModes * sizeof(REAL), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMalloc((void**) d_yModes, numAtoms * numModes * sizeof(REAL)));
		CUDA_CHECK(cudaMemcpy(d_yModes, protein->yModes(), numAtoms * numModes * sizeof(REAL), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMalloc((void**) d_zModes, numAtoms * numModes * sizeof(REAL)));
		CUDA_CHECK(cudaMemcpy(d_zModes, protein->zModes(), numAtoms * numModes * sizeof(REAL), cudaMemcpyHostToDevice));
	}

	typename DeviceProtein<REAL>::Desc deviceDesc;
	deviceDesc.numAtoms = numAtoms;
	deviceDesc.xPos = d_xPos;
	deviceDesc.yPos = d_yPos;
	deviceDesc.zPos = d_zPos;
	deviceDesc.type = d_type;
	deviceDesc.mappedType = d_mappedType;
	deviceDesc.charge = d_charge;
	deviceDesc.numModes = numModes;
	deviceDesc.xModes = d_xModes;
	deviceDesc.yModes = d_yModes;
	deviceDesc.zModes = d_zModes;

	std::shared_ptr<DeviceProtein<REAL>> deviceProtein = std::make_shared<DeviceProtein<REAL>>();
	deviceProtein->desc = deviceDesc;
	deviceProtein->hostResc = deviceDesc;

	return deviceProtein;
}

template<typename REAL>
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceProtein<REAL>> protein, deviceId_t deviceId) {
	ASSERT(deviceId >= 0);
	checkDeviceIdAndSetCurrent(deviceId);

	auto& hostResc = protein->hostResc;

	CUDA_CHECK(cudaFree(hostResc.xPos));
	CUDA_CHECK(cudaFree(hostResc.yPos));
	CUDA_CHECK(cudaFree(hostResc.zPos));
	CUDA_CHECK(cudaFree(hostResc.charge));
	CUDA_CHECK(cudaFree(hostResc.type));
	CUDA_CHECK(cudaFree(hostResc.mappedType));
	if (hostResc.numModes != 0) {
		CUDA_CHECK(cudaFree(hostResc.xModes));
		CUDA_CHECK(cudaFree(hostResc.yModes));
		CUDA_CHECK(cudaFree(hostResc.zModes));
	}
}

template<typename REAL>
std::shared_ptr<DeviceGridUnion<REAL>> DeviceDataConfigurator::attach(const std::shared_ptr<GridUnion<REAL>> grid, deviceId_t deviceId) {
	ASSERT(deviceId >= 0);
	checkDeviceIdAndSetCurrent(deviceId);
	std::shared_ptr<DeviceGridUnion<REAL>> d_grid = std::make_shared<DeviceGridUnion<REAL>>();

	d_grid->inner = attach(grid->inner);
	d_grid->outer = attach(grid->outer);
	d_grid->NL = attach(grid->NL);

	return d_grid;
}

template<typename REAL>
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceGridUnion<REAL>> d_grid, deviceId_t deviceId) {
	ASSERT(deviceId >= 0);
	checkDeviceIdAndSetCurrent(deviceId);

	detach(d_grid->inner);
	detach(d_grid->outer);
	detach(d_grid->NL);
}

template<typename REAL>
std::shared_ptr<DeviceIntrplGrid<REAL>> DeviceDataConfigurator::attach(const std::shared_ptr<IntrplGrid<REAL>> grid) {
	/** array of CUDA texture objects for built-in interpolation */
	cudaTextureObject_t* h_texArrayLin; /** located on host */
	cudaTextureObject_t* d_texArrayLin; /** located on device */
	/** array of CUDA texture objects for manual interpolation */
	cudaTextureObject_t* h_texArrayPt; 	/** located on host */
	cudaTextureObject_t* d_texArrayPt; 	/** located on device */
	cudaArray** h_cuArrayPtr;

	/** The following pointers are deleted when the grid gets detached from the device.
	 * They are needed to destroy device recourses (textures) and are kept in the hostResc object
	 * below.*/
	h_cuArrayPtr = new cudaArray*[grid->numTypes()];
	h_texArrayLin = new cudaTextureObject_t[grid->numTypes()];
	h_texArrayPt = new cudaTextureObject_t[grid->numTypes()];

	// Allocate CUDA array in device memory
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	auto dimN = grid->dimN();
	struct cudaExtent cuExtent = make_cudaExtent(dimN.x, dimN.y, dimN.z);
	// For cudaMalloc3DArray: width range in elements.
	// For cudaMalloc3D(...): width range in bytes.

	// Specify texture object parameters
	cudaTextureDesc texDescLin; // for built in device interpolation
	memset(&texDescLin, 0, sizeof(cudaTextureDesc));
	texDescLin.addressMode[0] = cudaAddressModeBorder; // return 0 if out of bounds
	texDescLin.addressMode[1] = cudaAddressModeBorder;
	texDescLin.addressMode[2] = cudaAddressModeBorder;
	texDescLin.filterMode = cudaFilterModeLinear;
	texDescLin.readMode = cudaReadModeElementType;
	texDescLin.normalizedCoords = false;

	cudaTextureDesc texDescPt; // for manual interpolation kernel
	memset(&texDescPt, 0, sizeof(cudaTextureDesc));
	texDescPt.addressMode[0] = cudaAddressModeBorder; // return 0 if out of bounds
	texDescPt.addressMode[1] = cudaAddressModeBorder;
	texDescPt.addressMode[2] = cudaAddressModeBorder;
	texDescPt.filterMode = cudaFilterModePoint;
	texDescPt.readMode = cudaReadModeElementType;
	texDescPt.normalizedCoords = false;

	for (unsigned i = 0; i<grid->numTypes(); i++) {
		cudaArray* &cuArray = h_cuArrayPtr[i];
		CUDA_CHECK(cudaMalloc3DArray(&cuArray, &channelDesc, cuExtent, cudaChannelFormatKindFloat));

		// copy data to 3D array
		cudaMemcpy3DParms copyParams;
		memset(&copyParams, 0, sizeof(copyParams));
		void* gridPtr = (void*)grid->getHostGridPtr(i);
		copyParams.srcPtr   = make_cudaPitchedPtr(gridPtr, cuExtent.width*sizeof(float4), cuExtent.width, cuExtent.height);
		copyParams.dstArray = cuArray;
		copyParams.extent   = cuExtent;
		copyParams.kind     = cudaMemcpyHostToDevice;
		CUDA_CHECK(cudaMemcpy3D(&copyParams));

		// Specify resource
		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = cuArray;

		// Create texture objects
		cudaTextureObject_t &texObjLin = h_texArrayLin[i];
		texObjLin = (long long)NULL;
		CUDA_CHECK(cudaCreateTextureObject(&texObjLin, &resDesc, &texDescLin, NULL));

		cudaTextureObject_t &texObjPt = h_texArrayPt[i];
		texObjPt = (long long)NULL;
		CUDA_CHECK(cudaCreateTextureObject(&texObjPt, &resDesc, &texDescPt, NULL));
	}

	CUDA_CHECK(cudaMalloc((void**)&d_texArrayLin, grid->numTypes()*sizeof(cudaTextureObject_t)));
	CUDA_CHECK(cudaMemcpy(d_texArrayLin, h_texArrayLin,grid->numTypes()*sizeof(cudaTextureObject_t), cudaMemcpyHostToDevice));

	CUDA_CHECK(cudaMalloc((void**)&d_texArrayPt, grid->numTypes()*sizeof(cudaTextureObject_t)));
	CUDA_CHECK(cudaMemcpy(d_texArrayPt, h_texArrayPt,grid->numTypes()*sizeof(cudaTextureObject_t), cudaMemcpyHostToDevice));

	typename DeviceIntrplGrid<REAL>::Desc desc;
	desc.dimN = grid->dimN();
	desc.dVox = grid->dVox();
	desc.dVox_inv = grid->dVox_inv();
	desc.voxelVol_inv = grid->voxelVol_inv();
 	desc.minDim = grid->pos();
 	desc.maxDim = grid->maxDim();
 	desc.texArrayLin = d_texArrayLin;
 	desc.texArrayPt = d_texArrayPt;

 	typename DeviceIntrplGrid<REAL>::HostResc hostResc;
 	hostResc.numArrays = grid->numTypes();
 	hostResc.cuArrayPtr = h_cuArrayPtr;
 	hostResc.h_texArrayLin = h_texArrayLin;
    hostResc.h_texArrayPt = h_texArrayPt;
    hostResc.d_texArrayLin = d_texArrayLin;
    hostResc.d_texArrayPt = d_texArrayPt;

    auto d_grid = std::make_shared<DeviceIntrplGrid<REAL>>();
    d_grid->desc = desc;
    d_grid->hostResc = hostResc;

    return d_grid;
}

template<typename REAL>
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceIntrplGrid<REAL>> grid) {
	auto& hostResc = grid->hostResc;
	/* Free interpolation grid resources */
	for (uint i = 0; i<hostResc.numArrays; i++) {
		CUDA_CHECK(cudaFreeArray(hostResc.cuArrayPtr[i]));
		CUDA_CHECK(cudaDestroyTextureObject(hostResc.h_texArrayLin[i]));
		CUDA_CHECK(cudaDestroyTextureObject(hostResc.h_texArrayPt[i]));
	}
	delete[] hostResc.cuArrayPtr;
	delete[] hostResc.h_texArrayLin;
	delete[] hostResc.h_texArrayPt;
	CUDA_CHECK(cudaFree(hostResc.d_texArrayLin));
	CUDA_CHECK(cudaFree(hostResc.d_texArrayPt));
}

template<typename REAL>
std::shared_ptr<DeviceNLGrid<REAL>> DeviceDataConfigurator::attach(const std::shared_ptr<NLGrid<REAL>> grid) {
	cudaChannelFormatDesc channelDesc =  cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindUnsigned);

	auto dimN = grid->dimN();
	struct cudaExtent cuExtent = make_cudaExtent(dimN.x, dimN.y, dimN.z);
	// For cudaMalloc3DArray: width range in elements.
	// For cudaMalloc3D(...): width range in bytes.

	// Specify texture object parameters
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder; // return 0 if out of bounds
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModePoint;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = false;

	cudaArray* cuArray;
	CUDA_CHECK(cudaMalloc3DArray(&cuArray, &channelDesc, cuExtent));

	// copy data to 3D array
	cudaMemcpy3DParms copyParams;
	memset(&copyParams, 0, sizeof(copyParams));
	copyParams.srcPtr = make_cudaPitchedPtr((void*) grid->grid(), cuExtent.width*sizeof(uint2), cuExtent.width, cuExtent.height);
	copyParams.dstArray = cuArray;
	copyParams.extent   = cuExtent;
	copyParams.kind     = cudaMemcpyHostToDevice;
	CUDA_CHECK(cudaMemcpy3D(&copyParams));

	// Specify resource
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = cuArray;

	// Create texture object
	cudaTextureObject_t texObj;
	texObj = (long long)NULL;
	CUDA_CHECK(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL));

	// Create device neighbor list
	unsigned* d_neighborList;
	CUDA_CHECK(cudaMalloc((void**)&d_neighborList, grid->neighborListSize()*sizeof(unsigned)));
	CUDA_CHECK(cudaMemcpy(d_neighborList, grid->neighborList(),grid->neighborListSize()*sizeof(unsigned), cudaMemcpyHostToDevice));

	typename DeviceNLGrid<REAL>::Desc desc;
	desc.dimN = grid->dimN();
	desc.dVox_inv  = grid->dVox_inv();
	desc.dPlateau2 = grid->dPlateau2();
	desc.minDim	= grid->minDim();
	desc.maxDim	= grid->maxDim();
	desc.tex = texObj;
	desc.neighborList = d_neighborList;

	typename DeviceNLGrid<REAL>::HostResc hostResc;
	hostResc.tex = texObj;
	hostResc.cuArray = cuArray;

	auto d_grid = std::make_shared<DeviceNLGrid<REAL>>();
    d_grid->desc = desc;
    d_grid->hostResc = hostResc;

	return d_grid;
}

template<typename REAL>
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceNLGrid<REAL>> grid) {
	auto& hostResc = grid->hostResc;
	/* Free NL grid resources */
	cudaVerify(cudaFreeArray(hostResc.cuArray));
	cudaVerify(cudaDestroyTextureObject(hostResc.tex));
}

template<typename REAL>
std::shared_ptr<DeviceParamTable<REAL>> DeviceDataConfigurator::attach(const std::shared_ptr<ParamTable<REAL>> paramTable, deviceId_t deviceId) {
	ASSERT(deviceId >= 0);
	checkDeviceIdAndSetCurrent(deviceId);

	const unsigned numTypes = paramTable->numTypes();

	typename ParamTable<REAL>::type_t* d_table;
	CUDA_CHECK(cudaMalloc((void**)&d_table, numTypes*numTypes*sizeof(typename ParamTable<REAL>::type_t)));
	CUDA_CHECK(cudaMemcpy(d_table, paramTable->table(), numTypes*numTypes*sizeof(typename ParamTable<REAL>::type_t), cudaMemcpyHostToDevice));


	typename DeviceParamTable<REAL>::Desc desc;
	desc.numTypes = numTypes;
	desc.shape = paramTable->potShape();
	desc.paramTable = d_table;

	typename DeviceParamTable<REAL>::HostResc hostResc;
	hostResc = desc;

	auto deviceParamTable = std::make_shared<DeviceParamTable<REAL>>();
	deviceParamTable->desc = desc;
	deviceParamTable->hostResc = hostResc;

	return deviceParamTable;
}

template<typename REAL>
void DeviceDataConfigurator::detach(const std::shared_ptr<DeviceParamTable<REAL>> d_table, deviceId_t deviceId) {
	ASSERT(deviceId >= 0);
	checkDeviceIdAndSetCurrent(deviceId);
	auto& hostResc = d_table->hostResc;

	CUDA_CHECK(cudaFree(hostResc.paramTable));
}


}  // namespace as

#endif // CUDA



#endif /* SRC_DEVICEDATACONFIGURATOR_TPP_ */
