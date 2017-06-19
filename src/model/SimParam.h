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

#ifndef SIMPARAM_H_
#define SIMPARAM_H_

#include "DataItem.h"

namespace as {
/*
 * The felec constant is the electrostatic energy, in kcal/mol,
 * between two electrons, at 1 A distance, in vacuum
 * The formula is:
 *  e**2 * NA * KC * Ang * kcal
 * where:
 *  e = charge of the electron, in Coulomb
 *  NA = Avogadro's number
 *  KC = Coulomb force constant
 *  Ang = the size of an Angstrom (10**-10 meter)
 *  kcal = the amount of Joules per kcal, 4184
 */
constexpr double FELEC = 332.053986;

enum class Dielec {
	constant,
	variable
};

template<typename REAL>
class SimParam : public DataItem {
	// Check if REAL is of floating-point type
	using real_t = typename std::enable_if<std::is_floating_point<REAL>::value, REAL>::type;
public:
	Dielec  dielec;		/** type of dielectric constant */
	real_t epsilon;		/** dielectric constant */
	real_t ffelec;		/** precomputed factor felec/epsilon */
};

}  // namespace





#endif /* SIMPARAM_H_ */
