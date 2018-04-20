/*
 * Types_6D.tpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_6D_MODES_TPP_
#define SRC_TYPES_6D_MODES_TPP_

#include <ostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "Types_6D_Modes.h"

namespace as {

template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_6D_Modes<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;


	outStream 	<< setw(w) << 0 << setw(w) << 0 << setw(w) << 0
				<< setw(w) << 0 << setw(w) << 0 << setw(w) << 0;
				for(int mode=0;mode<Common_Modes::numModesRec;mode++){
					outStream<< setw(w) << dof.modesRec[mode];
				}
	outStream	<< std::endl;

	outStream	<< setw(w) << dof._6D.ang.x << setw(w) << dof._6D.ang.y << setw(w) << dof._6D.ang.z;
	outStream	<< setw(w) << dof._6D.pos.x << setw(w) << dof._6D.pos.y << setw(w) << dof._6D.pos.z;
				for(int mode=0;mode<Common_Modes::numModesLig;mode++){
					outStream<< setw(w) << dof.modesLig[mode];
				}


	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result_6D_Modes<REAL> const& enGrad) {
	using namespace std;
	int precisionSetting = s.precision( );
	ios::fmtflags flagSettings = s.flags();

	s.setf(ios::scientific);
	s.precision(8);
	s << " Energy: " << enGrad._6D.E << endl;
	s.unsetf(ios::scientific);

	s.setf(ios::fixed);
	s.precision(3);
	s << setw(12) << enGrad._6D.E << endl;
	s.unsetf(ios::fixed);

	s.setf(ios::scientific);
	s.precision(8);
	int width = 20;


	s 	<< setw(width) << 0  << setw(width) << 0  << setw(width) << 0
		<< setw(width) << 0  << setw(width) << 0  << setw(width) << 0;
		for(int mode=0;mode<Common_Modes::numModesRec;mode++){
			s<< setw(width) << enGrad.modesRec[mode];
		}

	s 	<< setw(width) << enGrad._6D.ang.x  << setw(width) << enGrad._6D.ang.y  << setw(width) << enGrad._6D.ang.z
		<< setw(width) << enGrad._6D.pos.x  << setw(width) << enGrad._6D.pos.y  << setw(width) << enGrad._6D.pos.z;
		for(int mode=0;mode<Common_Modes::numModesLig;mode++){
			s<< setw(width) << enGrad.modesLig[mode];
		}

	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}

template<typename REAL>
void print_results( std::stringstream& os, DOF_6D_Modes<REAL> const & dof,  Result_6D_Modes<REAL> const& res){
	auto w = std::setw(20);
	os  << w << res._6D.E;
	os  << w <<  " " << w;
	os 	<< w << 0 << w << 0 << w << 0;
	os  << w << 0 << w << 0 << w << 0;
	for(int mode=0;mode<Common_Modes::numModesRec;mode++){
		os<< w << dof.modesRec[mode];
	}
	os  << w <<  " " << w;

	os	<< w << dof._6D.ang.x << w << dof._6D.ang.y << w << dof._6D.ang.z;
	os	<< w << dof._6D.pos.x << w << dof._6D.pos.y << w << dof._6D.pos.z;
	for(int mode=0;mode<Common_Modes::numModesLig;mode++){
		os<< w << dof.modesLig[mode];
	}
	os  << w <<  " " << w;
	os << std::endl;

}

template<typename REAL>
void print_results( std::stringstream& os,  Result_6D_Modes<REAL> const& res){
	auto w = std::setw(20);
	os  << w << res._6D.E;
	os  << w <<  " " << w;
	os 	<< w << 0 << w << 0 << w << 0;
	os  << w << 0 << w << 0 << w << 0;
	for(int mode=0;mode<Common_Modes::numModesRec;mode++){
		os<< w << res.modesRec[mode];
	}
	os  << w <<  " " << w;

	os	<< w << res._6D.ang.x << w << res._6D.ang.y << w << res._6D.ang.z;
	os	<< w << res._6D.pos.x << w << res._6D.pos.y << w << res._6D.pos.z;
	for(int mode=0;mode<Common_Modes::numModesLig;mode++){
		os<< w << res.modesLig[mode];
	}
	os  << w <<  " " << w;
	os << std::endl;

}


} // namespace as



#endif /* SRC_TYPES_6D_TPP_ */
