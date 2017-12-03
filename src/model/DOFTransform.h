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

#ifndef DOFTRANSFORM_H_
#define DOFTRANSFORM_H_

#include <vector>

#include "Vec3.h"
#include "RotMat.h"

namespace as {

template<typename REAL>
class DOF_6D;

template<typename REAL>
class DOF_6D_Modes;


	/*
	 ** @brief: This function performs an inplace transformation of the ligand coordinates assuming
	 ** that the receptor is always centered at the origin.
	 *
	 ** @input: [in] 	dof_rec: 	Vector of DOF-vectors of the receptor
	 ** 		[in+out]dof_lig: 	Vector of DOF-vectors of the ligand
	 ** 		[in]	pivot_rec:	Pivot of the receptor
	 ** 		[in]	pivot_lig:	Pivot of the ligand
	 ** 		[in]	center_rec: Are receptor coords centered?
	 ** 		[in]	center_lig: Are ligand coords centered?
	 **
	 ** 		pivots:	Vector of pivots of the receptor and the ligands
	 ** 				Can be obtained by calling asDB::readDOFHeader(...)
	 */
template<typename REAL>
void transformDOF_glob2rec(std::vector<DOF_6D<REAL>>& dof_rec, std::vector<DOF_6D<REAL>>& dof_lig,
			const Vec3<REAL>& pivot_rec, const Vec3<REAL>& pivot_lig,
			bool centered_rec, bool centered_lig);

template<typename REAL>
void transformDOF_glob2rec(std::vector<DOF_6D_Modes<REAL>>& dof_rec, std::vector<DOF_6D_Modes<REAL>>& dof_lig,
			const Vec3<REAL>& pivot_rec, const Vec3<REAL>& pivot_lig,
			bool centered_rec, bool centered_lig);

/* Does not work since torques cannot be rotated back -> force must be rotated before reduction */
//template<typename REAL, typename DOF, typename EnGrad>
//void transformEnGrad_rec2glob(const std::vector<DOF>& dof_rec, std::vector<EnGrad>& enGrad_lig);

} // namespace


#endif /* DOFTRANSFORM_H_ */
