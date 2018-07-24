/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2018 Glenn Glashagen
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

#include "monteCarlo.h"
#include <Eigen/Dense>
#include "RotMat.h"
#include "matrixFunctions.h"
using std::cerr;
using std::endl;

//as::MCSolver::Options as::MCSolver::settings;

extern "C" void minfor_(void* FortranSmuggler_ptr, int const& maxFunEval, int const& numModesRec, int const& numModesLig,
		double const* state);

namespace as {
/*
 ** @brief: creates a random rotation axis that is uniformly distributed on a spheres surface
 *
 ** @input: rand0 & rand1 are random numbers in the range [0, 1].
 *
 ** @output: axis[3], vector of length 1.0 .
 */
double PIx2 = 2*M_PI;
inline void randAxis(const double& rand0, const double& rand1, double axis[3]) {
	double phi = PIx2 * rand1;
	double theta = std::acos(2*rand0 - 1);
	double stheta = std::sin(theta);

	axis[0] = stheta*cos(phi);
	axis[1] = stheta*sin(phi);
	axis[2] = cos(theta);
}

/*
 ** @brief: creates a rotation matrix given an axis and an angle
 */
template<typename REAL>
inline void rotMatFromAxisAng(const double axis[3], const double ang, RotMat<REAL>& rotMat) {
	double c = std::cos(ang);
    double s = std::sin(ang);
    double t = 1.0 - c;

    rotMat[0] = c + axis[0]*axis[0]*t;
    rotMat[4] = c + axis[1]*axis[1]*t;
    rotMat[8] = c + axis[2]*axis[2]*t;


    double tmp1 = axis[0]*axis[1]*t;
    double tmp2 = axis[2]*s;
    rotMat[1] = tmp1 - tmp2;
    rotMat[3] = tmp1 + tmp2;

    tmp1 = axis[0]*axis[2]*t;
    tmp2 = axis[1]*s;
    rotMat[2] = tmp1 + tmp2;
    rotMat[6] = tmp1 - tmp2;


    tmp1 = axis[1]*axis[2]*t;
    tmp2 = axis[0]*s;
    rotMat[5] = tmp1 - tmp2;
    rotMat[7] = tmp1 + tmp2;
}


static double maxDist = 1;
static double maxAng = 3;
static double kT = 10;

/* Initialize random number generators */
static std::default_random_engine generator;
static std::uniform_real_distribution<double> distribution(0.0, 1.0);

template <typename REAL>
//inline void randomStep (const Vector& oldDOF, Vector& newDOF)
inline void randomStep (const double* oldDOF, double* newDOF)
{
//	int numdis = 20;
//	static int count = 0;
//	if (++count < numdis) {
//		std::cout << "#" << count << std::endl;
//	}
	/* get 6 random numbers: 3 pos. + 3 ang. */
	double r[6];
	for (unsigned k = 0; k < 6; ++k) {
		r[k] = distribution(generator);
	}

//	if (count < numdis) {
//		std::cout << "rand "<< r[0] << " "<< r[1] << " "<< r[2] << " "<< r[3] << " "<< r[4] << " "<< r[5] << std::endl;
//	}

	/********** Apply random rotation **********/

	/* create random axis */
	double rAxis[3];
	randAxis(r[0], r[1], rAxis);
//	if (count < numdis) {
//		std::cout << "rAxis "<< rAxis[0] << " "<< rAxis[1] << " "<< rAxis[2] << " "<< r[3] << " "<< r[4] << " "<< r[5] << std::endl;
//	}

	/* calculate rotation matrix according to axis and angle */
	double ang = (2.0*r[2]-1.0)*maxAng; // rad !!!
	RotMat<REAL> randRotMat;

	rotMatFromAxisAng(rAxis, ang, randRotMat);
	//if (count < numdis) {
	//	std::cout << "randRotMat " << randRotMat << endl;
	//}

	/* calculate current rotation matrix accoring to euler angles */

	double phi, ssi, rot;

//	phi = (double) oldDOF(0);
//	ssi = (double) oldDOF(1);
//	rot = (double) oldDOF(2);

	phi = (double) oldDOF[0];

	ssi = (double) oldDOF[1];
	rot = (double) oldDOF[2];

	RotMat<REAL> currRotMat = euler2rotmat<REAL>(phi, ssi, rot);
//	if (count < numdis) {
//		std::cout << "currRotMat " << currRotMat << endl;
//	}



	/* Multiply matrices and retain new angels */
	RotMat<REAL> newRotmat;
	newRotmat = randRotMat*currRotMat;
	rotmat2euler(newRotmat, phi, ssi, rot);
//	newDOF(0) = (float)phi;
//	newDOF(1) = (float)ssi;
//	newDOF(2) = (float)rot;

	newDOF[0] = (float)phi;
	newDOF[1] = (float)ssi;
	newDOF[2] = (float)rot;

	//std::cout << "phi" << (double)phi <<" "<< ssi << "  "<<rot << std::endl ;
	/********** Apply random displacement **********/

	double dx = (2.0*r[3] - 1.0)*maxDist;
	double dy = (2.0*r[4] - 1.0)*maxDist;
	double dz = (2.0*r[5] - 1.0)*maxDist;

//	cout << "#" << count << endl;
//	cout << rAxis[0] << " " << rAxis[1] << " " << rAxis[2] << " " << ang * 180.0 / M_PI << endl;
//	cout << dx << " " << dy << " " << dz << endl;
//	count++;

//	newDOF(3) = oldDOF(3) + dx;
//	newDOF(4) = oldDOF(4) + dy;
//	newDOF(5) = oldDOF(5) + dz;

	newDOF[0] = oldDOF[0] + dx;
	newDOF[1] = oldDOF[1] + dy;
	newDOF[2] = oldDOF[2] + dz;

//	if (count < numdis) {
//		std::cout << "old " << oldDOF << std::endl;
//		std::cout << "new " << newDOF << std::endl;
//	}
}


//void MC_accept(Vector& oldDOF, double& oldEnGrad, Vector &newDOF, double& newEnGrad) {
	void MC_accept(double* oldDOF, double& oldEnGrad, double *newDOF, double& newEnGrad) {
	double newEnergy = newEnGrad;
	double oldEnergy = oldEnGrad;

//	bool accepted = false;
	/* Metropolis Criterion */
	if (newEnergy <= oldEnergy) {
		oldEnGrad = newEnGrad;
		std::copy( newDOF, newDOF+6 , oldDOF);
		//oldDOF = newDOF;
//		accepted = true;
	} else {
		double r = distribution(generator);
		if (r < std::exp(-(newEnergy - oldEnergy)/kT)) {
			oldEnGrad = newEnGrad;
			std::copy( newDOF, newDOF+6 , oldDOF);
			//oldDOF = newDOF;
//			accepted = true;
		}
	}

//	int numdis = 20;
//	static int count = 0;
//	if (++count < numdis) {
//		std::cout << "accepted=" << accepted << " newE="<< newEnergy << " oldE=" << oldEnergy << std::endl;
//	}
}

void MCSolver::run(push_type& energyAndGradients) {
	/* Create Smuggler */
	//MCSolver::FortranSmuggler smuggler(ca, state, objective);

	/* create and fill state array */

	double energy = 0;
	ObjGrad objGrad_next;

	double* oldDOF = (double*) malloc(6 * sizeof(double));
	double* newDOF;

	newDOF = state.data();

	energyAndGradients();
	//std::cout << state << oldDOF << std::endl;
	/* the accepted values are stored in old variables !!! */
	//oldDof
	for ( int i = 0; i< 100; ++i){

		MC_accept(oldDOF, energy, newDOF, objective.obj);

		/* calulate new trial configuration */
		randomStep<double>(oldDOF, newDOF);
//		std::cout << "\nnew E: "<< energy << std::endl;
//		std::cout << "old DOF: "<< oldDOF[0] << oldDOF[1]<< oldDOF[2]<< oldDOF[3]<< oldDOF[4]<< oldDOF[5] << std::endl;
//		std::cout << "new DOF: "<< newDOF[0] << newDOF[1]<< newDOF[2]<< newDOF[3]<< newDOF[4]<< newDOF[5] << std::endl;



		energyAndGradients();
	}







	/* Algorithm 6.1, J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 140 */
//	Vector x_curr = getState();
//	Vector x_next;
//	const unsigned DIM = x_curr.rows();
//	Matrix H = Matrix::Identity(DIM, DIM);
//	ObjGrad objGrad_curr;
//	ObjGrad objGrad_next;
//
//
//	Vector p;
//	Vector s;
//	Vector y;
//
//	unsigned iter = 0;
//	double dObj = 0.0;
//
//
//	OBJGRAD(x_curr, objGrad_curr);
//	do {
//
//		if (stats) {
//			++statistic.numIter;
//		}
//
//		assert(std::isnan(objGrad_curr.obj) == false);
//		p = -1.0 * H * objGrad_curr.grad;
//
//		switch (search_type) {
//		case WolfeRule:
//			linesearch_WolfeRule(ca, x_curr, objGrad_curr, p, dObj, /*out*/x_next, /*out*/objGrad_next);
//			break;
//		default:
//			cerr << "Error: linesearch type unspecified" << endl;
//		}
//
//
//		s = x_next - x_curr;
//		y = objGrad_next.grad - objGrad_curr.grad;
//		dObj = objGrad_next.obj - objGrad_curr.obj;
//
//
//		double ys = y.dot(s);
//
//		if (fabs(dObj) < settings.dObjTol) {
//			if (stats) {
//				statistic.convergence = BFGSStatistic::finitePrec;
//			}
//			break;
//		}
//		if (s.norm() < settings.dxTol ) {
//			if (stats) {
//				statistic.convergence = BFGSStatistic::finitePrec;
//			}
//			break;
//		}
//
//		if (ys <= 0.0) {
//			/* damped BFGS-update taken from */
//			/* Al-Baali, Mehiddin, Lucio Grandinetti, and Ornella Pisacane.
//			 * "Damped techniques for the limited memory BFGS method for large-scale optimization."
//			 * Journal of Optimization Theory and Applications 161.2 (2014): 688-699. */
//			/* and */
//			/* Procedure 18.2, J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 537 */
//
//			Matrix B = H.inverse();
//
//			double theta;
//			Vector Bs = B*s;
//			double sBs = s.dot(Bs);
//			if (ys >= 0.2 * sBs) {
//				theta = 1.0; // skip updating
//			} else {
//				theta = 0.8 * sBs / (sBs - ys);
//			}
//			y = theta * y  +  (1.0-theta)*Bs;
//			assert(!std::isnan(theta) || !std::isinf(theta));
//
//			ys = y.dot(s);
//		}
//
//		assert(ys > 0.0);
//
//		if (ys < settings.illTol) {
//			if (stats) {
//				statistic.convergence = BFGSStatistic::illconditioned;
//			}
//			break;
//		}
//
//		const double rho = 1.0 / ys;
//
//		if (iter == 0) {
//			H = ((y.dot(s)) / (y.dot(y)) * Matrix::Identity(DIM, DIM));
//		}
//
//		H = (Matrix::Identity(DIM, DIM) - rho*s*y.transpose()) * H * (Matrix::Identity(DIM, DIM) - rho*y*s.transpose()) + rho*s * s.transpose();
//
//		x_curr = x_next;
//		objGrad_curr = objGrad_next;
//		iter++;
//
//	} while ((objGrad_next.grad.lpNorm<Eigen::Infinity>() > settings.gradTol) && (iter < settings.maxIter));
//
//	if (stats) {
//		statistic.gradNorm = objGrad_curr.grad.lpNorm<Eigen::Infinity>();
//		if (statistic.gradNorm <= settings.gradTol) {
//			statistic.convergence = BFGSStatistic::gradTolerance;
//		} else if (iter >= settings.maxIter) {
//			statistic.convergence = BFGSStatistic::maxIter;
//		}
//	}

	return;
}

} // namespace


// Call back function for fortran to access the class.  This is a free function
// that is not a member or friend of MyClass, so it can only access public
// member variables and functions of MyClass.  BaseClass::operator() is such a
// public member function.

