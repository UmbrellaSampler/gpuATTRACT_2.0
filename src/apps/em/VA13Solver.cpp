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
#include <iostream>

#include "VA13Solver.h"

using std::cerr;
using std::endl;

as::VA13Solver::Options as::VA13Solver::settings;

extern "C" void minfor_(void* FortranSmuggler_ptr, int const& maxFunEval, int const& minTrans,int const& minRot,int const& minMode,int const& numModesRec, int const& numModesLig,
		double const* state);

namespace as {


void VA13Solver::run(push_type& ca) {
	/* Create Smuggler */
	VA13Solver::FortranSmuggler smuggler(ca, state, objective, trackedStates, trackedGrads, settings);
	/* create and fill state array */

	double state_array[state.rows()];
	for (int i = 0; i < state.rows(); ++i) {
		state_array[i] = state(i);

	}

	minfor_(&smuggler, settings.maxFunEval,
			settings.minimizeTranslation,
			settings.minimizeRotation,
			settings.minimizeModes,
			Common_Modes::numModesRec,
			Common_Modes::numModesLig, state_array);
}


} // namespace


// Call back function for fortran to access the class.  This is a free function
// that is not a member or friend of MyClass, so it can only access public
// member variables and functions of MyClass.  BaseClass::operator() is such a
// public member function.
extern "C" void energy_for_fortran_to_call_(void* FortranSmuggler_ptr, double state_ptr[], double* energy, double grad[])
{
   // Cast to BaseClass.  If the pointer isn't a pointer to an object
   // derived from BaseClass, then the world will end.
	as::VA13Solver::FortranSmuggler* smuggler = static_cast<as::VA13Solver::FortranSmuggler*>(FortranSmuggler_ptr);

	/* set the state */
//	std::cout << "state ";
	as::Vector& state = smuggler->state_ref();
	for (int i = 0; i < state.rows(); ++i) {
		state(i) = state_ptr[i];
	//	std::cout <<std::setprecision(10) <<" "<<state(i);
	}
//	std::cout << endl;
	/* call coroutine to break execution here until energy and gradients are available */
	smuggler->call_coro();

	/* get the state */
	as::ObjGrad& objGrad = smuggler->objective_ref();
	*energy = objGrad.obj;
	//std::cout << "grad "<< objGrad.obj<<  " ";
	for (int i = 0; i < objGrad.grad.rows(); ++i) {
		grad[i] = objGrad.grad(i);
	}


}


extern "C" void state_tracker_(void* FortranSmuggler_ptr, double state_ptr[], double* energy, double grad[], int * size)
{

	as::VA13Solver::FortranSmuggler* smuggler = static_cast<as::VA13Solver::FortranSmuggler*>(FortranSmuggler_ptr);

	as::VA13Solver::Options settings = smuggler->getOptions();
	if (settings.trackGradients){
		std::vector<float> grads;
		grads.push_back(*energy);
		for (int i = 0; i < *size; ++i)
		{
			grads.push_back(grad[i]);
		}
		smuggler->push_grad(grads);
	}
	if (settings.trackStates){
		std::vector<float> state;
		for (int i = 0; i < *size; ++i)
		{
			state.push_back(state_ptr[i]);
		}
		smuggler->push_state(state);
	}
}
