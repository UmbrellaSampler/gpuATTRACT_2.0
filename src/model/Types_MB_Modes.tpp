/*
 * Types_MB.tpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_MB_MODES_TPP_
#define SRC_TYPES_MB_MODES_TPP_

#include <ostream>
#include <iomanip>
#include <sstream>
#include "Types_MB_Modes.h"

namespace as {


template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_Modes<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;


	outStream 	<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z
				<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z;
				for(int mode=0;mode<MODES_MAX_NUMBER;mode++){
					outStream<< setw(w) << dof.modes[mode];
				}


	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_MB_Modes<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;

	for ( int i = 0; i < NUM_MAX_PROTEIN; ++i){
		outStream << dof.protein[i]<< std::endl;
	}



	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result_MB_Modes<REAL> const& enGrad) {
	using namespace std;
	int precisionSetting = s.precision( );
	ios::fmtflags flagSettings = s.flags();

	s.setf(ios::scientific);
	s.precision(8);
	s << " Energy: " << enGrad.E << endl;
	s.unsetf(ios::scientific);

	s.setf(ios::fixed);
	s.precision(3);
	s << setw(12) << enGrad.E << endl;
	s.unsetf(ios::fixed);

	s.setf(ios::scientific);
	s.precision(8);
	int width = 20;


	for ( int i = 0; i < NUM_MAX_PROTEIN; ++i){
			s << enGrad.protein[i]<< std::endl;
		}

	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}


} // namespace as



#endif /* SRC_TYPES_MB_TPP_ */
