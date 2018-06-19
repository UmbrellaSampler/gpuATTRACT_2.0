/*
 * DeviceProtein.h
 *
 *  Created on: May 29, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEPROTEIN_H_
#define SRC_DEVICEPROTEIN_H_

#ifdef CUDA
#include "DeviceItem.h"
#include "Types_6D_Config.h"
/**
 * DeviceProtein represents a Protein on a device.
 */
namespace as {

template <typename REAL>
class DeviceProtein : public DeviceItem {
public:

	struct Desc {
		REAL *xPos = nullptr;	/** Cartesian coordinates in cm-frame*/
		REAL *yPos = nullptr;
		REAL *zPos = nullptr;
		REAL* xModes = nullptr; /** normal mode deformation vectors */
		REAL* yModes = nullptr;
		REAL* zModes = nullptr;
		REAL* charge;	/** charge of the atoms/particle */
		unsigned* type = nullptr; 	/** atom type */
		unsigned* mappedType = nullptr;
		unsigned numMappedTypes;
		unsigned numModes; /** number of modes */
		unsigned numAtoms; 	/** number of atoms/particles */
		REAL modeForce[MODES_MAX_NUMBER];
		bool isOrigin;
	};

	using HostResc = Desc;

	Desc desc;
	HostResc hostResc;
};

template<typename REAL>
using d_Protein = typename DeviceProtein<REAL>::Desc;

}

#endif // CUDA

#endif /* SRC_DEVICEPROTEIN_H_ */
