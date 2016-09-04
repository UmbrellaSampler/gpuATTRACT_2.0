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

#ifndef SRC_PROTEIN_TPP_
#define SRC_PROTEIN_TPP_

#include <stdexcept>

#include "Protein.h"

//todo: remove
#include <iostream>
#include <iomanip>

namespace as {

template<typename REAL>
Protein<REAL>::Protein() :
	_tag(),
	_numAtoms(0), _pivot(0,0,0),
	_pos(nullptr),
	_type(nullptr),
	_numMappedTypes(0),
	_mappedTypes(nullptr),
	_charge(nullptr),
	_numModes(0), _modes(nullptr) {};

template<typename REAL>
Protein<REAL>::~Protein() {
	delete[] _pos;
	delete[] _charge;
	delete[] _type;
	delete[] _mappedTypes;
	delete[] _modes;
}

template<typename REAL>
REAL* Protein<REAL>::getOrCreatePosPtr() {
	if (_pos == nullptr) {
		if (_numAtoms == 0) {
			throw std::runtime_error("getOrCreatePosPtr(): the number of atoms must be set before");
		}
		_pos = new REAL[3*_numAtoms];
	}
	return _pos;
}

template<typename REAL>
unsigned* Protein<REAL>::getOrCreateTypePtr() {
	if (_type == nullptr) {
		if (_numAtoms == 0) {
			throw std::runtime_error("getOrCreateTypePtr(): the number of atoms must be set before");
		}
		_type = new type_t[_numAtoms];
	}
	return _type;
}

template<typename REAL>
auto Protein<REAL>::getOrCreateMappedPtr() -> type_t* {
	if (_mappedTypes == nullptr) {
		if (_numAtoms == 0) {
			throw std::runtime_error("getOrCreateMappedPtr(): the number of atoms must be set before");
		}
		if (_numMappedTypes == 0) {
			throw std::runtime_error("getOrCreateMappedPtr(): the number of mapped types must be set before");
		}
		_mappedTypes = new type_t[_numAtoms*_numMappedTypes];
	}
	return _mappedTypes;
}

template<typename REAL>
REAL* Protein<REAL>::getOrCreateChargePtr() {
	if (_charge == nullptr) {
		if (_numAtoms == 0) {
			throw std::runtime_error("getOrCreateChargePtr(): the number of atoms must be set before");
		}
		_charge = new REAL[_numAtoms];
	}
	return _charge;
}

template<typename REAL>
REAL* Protein<REAL>::getOrCreateModePtr() {
	if (_modes == nullptr) {
		if (_numAtoms == 0) {
			throw std::runtime_error("getOrCreateModePtr(): the number of atoms must be set before");
		}

		if (_numModes == 0) {
			throw std::runtime_error("getOrCreateModePtr(): the number of modes must be set before");
		}
		_modes = new REAL[3*_numAtoms*_numModes];
	}
	return _modes;
}

template<typename REAL>
void Protein<REAL>::pivotize(vec3_t pivot) {
	if (_pivot != vec3_t(0,0,0)) {
		undoPivoting();
	}
	_pivot = pivot;
	if (_pivot != vec3_t(0,0,0)) {
		for (unsigned i = 0; i < _numAtoms; ++i) {
			xPos()[i] -= _pivot.x;
			yPos()[i] -= _pivot.y;
			zPos()[i] -= _pivot.z;
		}
	}
}

template<typename REAL>
void Protein<REAL>::auto_pivotize() {
	if (_pivot != vec3_t(0,0,0)) {
		undoPivoting();
	}
	vec3_t pivot(0,0,0);
	for (unsigned i = 0; i < _numAtoms; ++i) {
		pivot.x += xPos()[i];
		pivot.y += yPos()[i];
		pivot.z += zPos()[i];
	}
	pivot /= static_cast<double>(_numAtoms);
	pivotize(pivot);
}

template<typename REAL>
void Protein<REAL>::undoPivoting() {
	for (unsigned i = 0; i < _numAtoms; ++i) {
		xPos()[i] += _pivot.x;
		yPos()[i] += _pivot.y;
		zPos()[i] += _pivot.z;
	}
	_pivot = vec3_t(0,0,0);
}

//TODO: move print to free functions
template<typename REAL>
void Protein<REAL>::print(int numEL) const {
	using namespace std;
	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();
	cout.setf(ios::dec | ios::showpoint | ios::showpos);
	cout.precision(6);

	int w = 13;
//	outStream 	<< setw(w) << "DOF"
//				<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z
//				<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z;

	cout << setw(5) << "#" << setw(w) << "X" << setw(w) << "Y" << setw(w) << "Z" << setw(6) << "TYPE" << setw(w) << "CHARGE" << endl;
	for(unsigned i = 0; i < numEL; ++i) {
		cout << setw(5) << i+1
			 << setw(w) << xPos()[i]
		     << setw(w) << yPos()[i]
		     << setw(w) << zPos()[i]
		     << setw(6) << type()[i]
		     << setw(6) << mappedType()[i]
			 << setw(w) << charge()[i]
			 << endl;
	}

	cout.precision(precisionSetting);
	cout.flags(flagSettings);
}

} // namespace as

#endif






