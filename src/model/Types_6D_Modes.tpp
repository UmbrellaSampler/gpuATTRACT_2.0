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
std::ostream& operator<<(std::ostream& outStream, DOF2_6D_Modes<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;


	outStream 	<< setw(w) << "DOF" << setw(w) << dof.rec.pos.x << setw(w) << dof.rec.pos.y << setw(w) << dof.rec.pos.z
				<< setw(w) << dof.rec.ang.x << setw(w) << dof.rec.ang.y << setw(w) << dof.rec.ang.z;
				for(int mode=0;mode<dof.rec.numModes;mode++){outStream<< setw(w) << dof.rec.modes[mode];}

	outStream 	<< setw(w) << "DOF" << setw(w) << dof.lig.pos.x << setw(w) << dof.lig.pos.y << setw(w) << dof.lig.pos.z
				<< setw(w) << dof.lig.ang.x << setw(w) << dof.lig.ang.y << setw(w) << dof.lig.ang.z;
				for(int mode=0;mode<dof.lig.numModes;mode++){outStream<< setw(w) << dof.lig.modes[mode];}

	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result2_6D_Modes<REAL> const& enGrad) {
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


	s << " Gradients: " << setw(width) << enGrad.rec.ang.x  << setw(width) << enGrad.rec.ang.y  << setw(width) << enGrad.rec.ang.z
				<< setw(width) << enGrad.rec.pos.x  << setw(width) << enGrad.rec.pos.y  << setw(width) << enGrad.rec.pos.z;
				for(int mode=0;mode<enGrad.rec.numModes;mode++){s<< setw(width) << enGrad.rec.modes[mode];}

	s << " Gradients: " << setw(width) << enGrad.lig.ang.x  << setw(width) << enGrad.lig.ang.y  << setw(width) << enGrad.lig.ang.z
				<< setw(width) << enGrad.lig.pos.x  << setw(width) << enGrad.lig.pos.y  << setw(width) << enGrad.lig.pos.z;
				for(int mode=0;mode<enGrad.lig.numModes;mode++){s<< setw(width) << enGrad.lig.modes[mode];}


	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}

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


} // namespace as



#endif /* SRC_TYPES_6D_TPP_ */
