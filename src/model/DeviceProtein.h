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
#include "Vec3.h"
/**
 * DeviceProtein represents a Protein on a device.
 */
namespace as {

template <typename REAL>
class DeviceProtein : public DeviceItem {
public:

	struct Desc {
		unsigned numAtoms; 	/** number of atoms/particles */

		REAL *xPos = nullptr;	/** Cartesian coordinates in cm-frame*/
		REAL *yPos = nullptr;
		REAL *zPos = nullptr;
		REAL* charge;	/** charge of the atoms/particle */
		REAL* xModes = nullptr; /** normal mode deformation vectors */
		REAL* yModes = nullptr;
		REAL* zModes = nullptr;
		unsigned* type = nullptr; 	/** atom type */
		unsigned* mappedType = nullptr;
		REAL modeForce[MODES_MAX_NUMBER];
		Vec3<REAL> pivot;
		unsigned numModes; /** number of modes */
		unsigned numMappedTypes;
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
