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


	outStream 	<< setw(w) << "DOF"
				<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z
				<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z;
				for(int mode=0;mode<dof.numModes;mode++){outStream<< setw(w) << dof.modes[mode];}

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
	s << " Energy: " << enGrad.E << endl;
	s.unsetf(ios::scientific);

	s.setf(ios::fixed);
	s.precision(3);
	s << setw(12) << enGrad.E << endl;
	s.unsetf(ios::fixed);

	s.setf(ios::scientific);
	s.precision(8);
	int width = 20;


	s << " Gradients: "
			<< setw(width) << enGrad.ang.x  << setw(width) << enGrad.ang.y  << setw(width) << enGrad.ang.z
			<< setw(width) << enGrad.pos.x  << setw(width) << enGrad.pos.y  << setw(width) << enGrad.pos.z;
	for(int mode=0;mode<enGrad.numModes;mode++){s<< setw(width) << enGrad.modes[mode];}
	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}

} // namespace as



#endif /* SRC_TYPES_6D_TPP_ */
