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

#ifndef META_H_
#define META_H_

#include <Eigen/Core>
#include "Types_6D.h"


#define OBJGRAD(dof, energy)	\
	do { 						\
		/*std::cout << "\t" << "request=" << Vector2extDOF(dof) << std::endl;*/ \
		state = dof;			\
		ca();					\
		energy = objective;		\
		/*std::cout << "\t" << "result=" << ObjGrad2extEnGrad(energy) << std::endl;*/						\
	} while(0)


#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x,y) ((x < y) ? y : x)
#endif

namespace as {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Scalar = Eigen::VectorXd::Scalar;

struct ObjGrad {
	double obj; // function value
	Vector grad; // gradients
};


template<typename FIRST, typename SECOND>
class TypesConverter {
	static FIRST toFirst(SECOND const&);
	static SECOND toSecond(FIRST const&);
};


template<typename REAL>
class TypesConverter<DOF_6D_Impl<REAL>, Vector> {
public:
	static DOF_6D_Impl<REAL> toFirst(Vector const& vec) {
		DOF_6D_Impl<REAL> dof;
		dof.ang.x = vec(0);
		dof.ang.y = vec(1);
		dof.ang.z = vec(2);
		dof.pos.x = vec(3);
		dof.pos.y = vec(4);
		dof.pos.z = vec(5);
		return dof;
	}

	static Vector toSecond(DOF_6D_Impl<REAL> const& dof) {
		Vector vec(6);
		vec  << dof.ang.x, dof.ang.y, dof.ang.z,
				dof.pos.x, dof.pos.y , dof.pos.z;
		return vec;
	}
};

template<typename REAL>
class TypesConverter<Result_6D_Impl<REAL>, ObjGrad> {
public:
	static Result_6D_Impl<REAL> toFirst(ObjGrad const& objGrad) {
		Result_6D_Impl<REAL> enGrad;
		enGrad.E_VdW = objGrad.obj;
		enGrad.E_El = 0.0; // TODO This is not true, we just lost the information
		enGrad.ang.x = objGrad.grad(0);
		enGrad.ang.y = objGrad.grad(1);
		enGrad.ang.z = objGrad.grad(2);
		enGrad.pos.x = objGrad.grad(3);
		enGrad.pos.y = objGrad.grad(4);
		enGrad.pos.z = objGrad.grad(5);
		return enGrad;
	}

	static ObjGrad toSecond (Result_6D_Impl<REAL> const& enGrad) {
		ObjGrad objGrad;
		objGrad.obj = enGrad.E_VdW + enGrad.E_El;
		objGrad.grad = Vector(6);
		// for ATTRACT multiply gradients by -1.0
		objGrad.grad  << -enGrad.ang.x, -enGrad.ang.y,  -enGrad.ang.z,
						 -enGrad.pos.x, -enGrad.pos.y , -enGrad.pos.z;
		return objGrad;
	}
};

} //namespace

#endif /* META_H_ */
