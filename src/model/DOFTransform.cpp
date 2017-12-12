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

#include <vector>
#include <cmath>
#include <cassert>

#include "nativeTypesMath.h"
#include "nativeTypesWrapper.h"
#include "DOFTransform.h"
#include "Types_6D.h"
#include "Types_6D_Modes.h"
#include "matrixFunctions.h"
#include <iostream>

namespace as {


template<typename REAL>
void transformDOF_glob2rec(DOF_6D<REAL>& dof_rec, DOF_6D<REAL>& dof_lig,
			const Vec3<REAL>& pivot_rec, const Vec3<REAL>& pivot_lig,
			bool const& centered_rec, bool const& centered_lig) {
	using namespace std;
	using vec3_t = Vec3<REAL>;

	// center ligand dofs by the respective ligand pivot
	if (!centered_lig) {
		dof_lig.pos.x += pivot_lig.x;
		dof_lig.pos.y += pivot_lig.y;
		dof_lig.pos.z += pivot_lig.z;
	}

	/* shift ligand dofs by receptor pivot */
	if (!centered_rec) {
		if (!(pivot_rec == Vec3<REAL> (0.0))) {
			dof_lig.pos.x -= pivot_rec.x;
			dof_lig.pos.y -= pivot_rec.y;
			dof_lig.pos.z -= pivot_rec.z;
		}
	}



	/* rotate ligand into receptor frame and shift ligand by receptor dofs*/
	const vec3_t& pos_rec = dof_rec.pos;

	if (pos_rec.x != 0.0 || pos_rec.y != 0.0 || pos_rec.z != 0.0 ) {
		vec3_t& pos_lig = dof_lig.pos;
		pos_lig = pos_lig - pos_rec;
	}
	const vec3_t& ang_rec = dof_rec.ang;
	if (ang_rec.x != 0.0 || ang_rec.y != 0.0 || ang_rec.z != 0.0 ) {

		vec3_t& ang_lig = dof_lig.ang;
		vec3_t& pos_lig = dof_lig.pos;

		RotMat<REAL> mat_rec_inv = euler2rotmat(-ang_rec.x, -ang_rec.y, -ang_rec.z);
		RotMat<REAL> mat_lig = euler2rotmat(ang_lig.x, ang_lig.y, ang_lig.z);

		pos_lig = mat_rec_inv * pos_lig;
		RotMat<REAL> mat =  mat_rec_inv *mat_lig ;
		rotmat2euler(mat, ang_lig.x, ang_lig.y, ang_lig.z);
	}

	/* the receptor dofs can now considered to be zero */
	dof_rec.pos = vec3_t(0.0, 0.0, 0.0);
	dof_rec.ang = vec3_t(0.0, 0.0, 0.0);
}

template<typename REAL>
void transformDOF_glob2rec(std::vector<DOF_6D<REAL>>& dofs_rec, std::vector<DOF_6D<REAL>>& dofs_lig,
			const Vec3<REAL>& pivot_rec, const Vec3<REAL>& pivot_lig,
			bool centered_rec, bool centered_lig) {
	using namespace std;

	assert(dofs_rec.size() == dofs_lig.size());

	for (size_t i = 0; i < dofs_rec.size(); ++i) {
		transformDOF_glob2rec(dofs_rec[i], dofs_lig[i],
				pivot_rec, pivot_lig,
				centered_rec, centered_lig);
	}
}

template<typename REAL>
void transformDOF_glob2rec(std::vector<DOF_6D_Modes<REAL>>& dofs_rec, std::vector<DOF_6D_Modes<REAL>>& dofs_lig,
			const Vec3<REAL>& pivot_rec, const Vec3<REAL>& pivot_lig,
			bool centered_rec, bool centered_lig) {
	using namespace std;

	assert(dofs_rec.size() == dofs_lig.size());

	for (size_t i = 0; i < dofs_rec.size(); ++i) {
		transformDOF_glob2rec(dofs_rec[i]._6D, dofs_lig[i]._6D,
				pivot_rec, pivot_lig,
				centered_rec, centered_lig);
	}
}

template
void transformDOF_glob2rec(std::vector<DOF_6D<float>>& dof_rec, std::vector<DOF_6D<float>>& dof_lig,
			const Vec3<float>& pivot_rec, const Vec3<float>& pivot_lig,
			bool centered_rec, bool centered_lig);

template
void transformDOF_glob2rec(std::vector<DOF_6D<double>>& dof_rec, std::vector<DOF_6D<double>>& dof_lig,
			const Vec3<double>& pivot_rec, const Vec3<double>& pivot_lig,
			bool centered_rec, bool centered_lig);

template
void transformDOF_glob2rec(std::vector<DOF_6D_Modes<float>>& dof_rec, std::vector<DOF_6D_Modes<float>>& dof_lig,
			const Vec3<float>& pivot_rec, const Vec3<float>& pivot_lig,
			bool centered_rec, bool centered_lig);

template
void transformDOF_glob2rec(std::vector<DOF_6D_Modes<double>>& dof_rec, std::vector<DOF_6D_Modes<double>>& dof_lig,
			const Vec3<double>& pivot_rec, const Vec3<double>& pivot_lig,
			bool centered_rec, bool centered_lig);

//template<typename REAL, typename DOF, typename EnGrad>
//void transformEnGrad_rec2glob(const std::vector<DOF>& dof_rec, std::vector<EnGrad>& enGrad_lig) {
//	using namespace std;
//	using real3_t = typename TypeWrapper<REAL>::real3_t;
//
//	assert(dof_rec.size() == enGrad_lig.size());
//
//	/* rotate forces and torques in global frame */
//	for (unsigned j = 0; j < dof_rec.size(); ++j) {
//		const real3_t& ang_rec = dof_rec[j].ang;
//		if (ang_rec.x != 0.0 || ang_rec.y != 0.0 || ang_rec.z != 0.0 ) {
//			real3_t& force_lig = enGrad_lig[j].pos;
//			Vec3<REAL> force_lig_v(force_lig.x, force_lig.y, force_lig.z);
//
//			/* rotate forces */
//
//			RotMat<REAL> mat_rec = euler2rotmat(ang_rec.x, ang_rec.y, ang_rec.z);
//
//			force_lig_v = mat_rec * force_lig_v;
//			force_lig = make_real(force_lig_v[0], force_lig_v[1], force_lig_v[2]);
//
//			/* rotate torques */
//
//
//		}
//	}
//}

}  // namespace as


