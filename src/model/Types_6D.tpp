/*
 * Types_6D.tpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_6D_TPP_
#define SRC_TYPES_6D_TPP_

#include <ostream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include "Types_6D.h"

namespace as {


template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_6D<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;
	//outStream 	<< setw(w) << "DOF"
	outStream 	<< "#" << endl;
	outStream	<< setw(w) << "0" << setw(w) << "0"<< setw(w) << "0"<< setw(w) << "0"<< setw(w) << "0"<< setw(w) << "0"<<endl;
	outStream	<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z;
	outStream	<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z;
	outStream	<< endl;
	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result_6D<REAL> const& enGrad) {
	using namespace std;
	int precisionSetting = s.precision( );
	ios::fmtflags flagSettings = s.flags();

	s.setf(ios::scientific);
	s.precision(8);
	s << " Energy: " << enGrad.E << endl;
	s.unsetf(ios::scientific);

	s.setf(ios::fixed);
	s.precision(3);
	s << setw(12) << enGrad.E <<endl;
	//s << setw(12) << " ; ";// << endl;
	s.unsetf(ios::fixed);

	s.setf(ios::scientific);
	s.precision(8);
	int width = 20;
	s << " Gradients: ";
	s		<< setw(width) << enGrad.ang.x  << setw(width) << enGrad.ang.y  << setw(width) << enGrad.ang.z
			<< setw(width) << enGrad.pos.x  << setw(width) << enGrad.pos.y  << setw(width) << enGrad.pos.z;
	//s		<< " ; ";
	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}

template<typename REAL>
void print_results( std::stringstream& os,  DOF_6D<REAL> const & dof, Result_6D<REAL> const & res){
	auto w = std::setw(20);
	os  << w << res.E;
	os  << w <<  " " << w;
	os	<< w << dof.ang.x << w << dof.ang.y << w << dof.ang.z;
	os	<< w << dof.pos.x << w << dof.pos.y << w << dof.pos.z;
	os  << w <<  " " << w;
	os << std::endl;
}

template<typename REAL>
void print_results( std::stringstream& os,  Result_6D<REAL> const & res){
	auto w = std::setw(20);
	os  << w << res.E;
	os  << w <<  " " << w;
	os	<< w << res.ang.x << w << res.ang.y << w << res.ang.z;
	os	<< w << res.pos.x << w << res.pos.y << w << res.pos.z;
	os  << w <<  " " << w;
	os << std::endl;
}


} // namespace as



#endif /* SRC_TYPES_6D_TPP_ */
