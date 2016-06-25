/*
 * DeviceProtein.h
 *
 *  Created on: May 29, 2016
 *      Author: uwe
 */

#ifndef SRC_DEVICEPROTEIN_H_
#define SRC_DEVICEPROTEIN_H_

/**
 * DeviceProtein represents a Protein on a device.
 */
namespace as {

template <typename REAL>
class DeviceProtein;


template <typename REAL>
class DeviceProtein {
public:
	template <typename _REAL>
	struct Desc {
		unsigned numAtoms; 	/** number of atoms/particles */

		_REAL *xPos;	/** Cartesian coordinates in cm-frame*/
		_REAL *yPos;
		_REAL *zPos;

		unsigned* type; 	/** atom type */
		unsigned numMappedTypes;
		unsigned* mappedType;

		_REAL* charge;	/** charge of the atoms/particle */

		unsigned numModes; /** number of modes */
		_REAL* xModes; /** normal mode deformation vectors */
		_REAL* yModes;
		_REAL* zModes;
	};

	using HostResc = Desc<REAL>;

	Desc<REAL> desc;
	HostResc hostResc;
};

//template <typename REAL>
//template <typename _REAL>
//struct DeviceProtein<REAL>::Desc<_REAL> {
//	unsigned numAtoms; 	/** number of atoms/particles */
//
//	_REAL *xPos;	/** Cartesian coordinates in cm-frame*/
//	_REAL *yPos;
//	_REAL *zPos;
//
//	unsigned* type; 	/** atom type */
//	unsigned numMappedTypes;
//	unsigned* mappedType;
//
//	_REAL* charge;	/** charge of the atoms/particle */
//
//	unsigned numModes; /** number of modes */
//	_REAL* xModes; /** normal mode deformation vectors */
//	_REAL* yModes;
//	_REAL* zModes;
//};


}
#endif /* SRC_DEVICEPROTEIN_H_ */
