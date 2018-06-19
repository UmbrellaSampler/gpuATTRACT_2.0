/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include <string>
#include <map>
#include <algorithm>
#include <cassert>
#include <type_traits>

#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "TypeMap.h"
#include "DataItem.h"


namespace as {

template<typename REAL>
class Protein : public DataItem {
	// Check if REAL is of floating-point type
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;

public:
	using type_t = unsigned;
	using tag_t = std::string;

	Protein();
	virtual ~Protein();

	unsigned numAtoms() const{
		return _numAtoms;
	}

	REAL* xPos() const{
		return _pos;
	}

	REAL* yPos() const{
		return _pos + _numAtoms;
	}

	REAL* zPos() const{
		return _pos + 2*_numAtoms;
	}

	REAL* charge() const{
		return _charge;
	}

	type_t* type() const{
		return _type;
	}

	size_t numMappedTypes() {
		return _numMappedTypes;
	}

	type_t* mappedType(unsigned gridId = 0) const {
		assert(gridId < _numMappedTypes);
		return _mappedTypes + gridId*_numAtoms;
	}

	unsigned numModes() const{
		return _numModes;
	}

	REAL* xModes() const {
		return _modes;
	}

	REAL* yModes() const{
		return _modes + _numModes*_numAtoms;
	}

	REAL* zModes() const{
		return _modes + 2*_numModes*_numAtoms;
	}

	REAL* modeForce(){
		return _modeForceConstant;
	}

	tag_t tag() const{
		return _tag;
	}

	vec3_t pivot() {
		return _pivot;
	}

	REAL* getOrCreatePosPtr();

	type_t* getOrCreateTypePtr();

	type_t* getOrCreateMappedPtr();

	REAL* getOrCreateChargePtr();

	REAL* getOrCreateModePtr();

	REAL* getOrCreateModeForcePtr();

	void setTag(tag_t tag) {
		_tag = tag;
	}

	void setNumAtoms(unsigned num) {
		_numAtoms = num;
	}

	void setNumModes(unsigned num) {
		_numModes = num;
	}

	void setNumMappedTypes(unsigned num) {
		_numMappedTypes = num;
	}

	//TODO: refactor the pivotize functions to non-members
	void pivotize(vec3_t pivot);

	void auto_pivotize();

	void undoPivoting();

	void print(int numEl) const;

	void setOrigin( bool isOrigin ){
		_isOrigin = isOrigin;
	}

	bool getOrigin(){
		return _isOrigin;
	}

protected:
	std::string _tag;	/** identifier: filename (default) */
	unsigned _numAtoms; /** number of atoms/particles */


	vec3_t _pivot;	/** rotation pivot */
	REAL *_pos;	/** Cartesian coordinates in cm-frame (pivotized) */

	type_t* _type; 	/** atom type */
	unsigned _numMappedTypes;
	type_t* _mappedTypes; /* for receptor grid mapping */

	REAL* _charge;	/** charge of the atoms/particle */
	bool _isOrigin;

	unsigned _numModes; /** number of modes */
	REAL* _modes; /** normal mode deformation vectors */

	REAL* _modeForceConstant;
};

}


#endif /* PROTEIN_H_ */

