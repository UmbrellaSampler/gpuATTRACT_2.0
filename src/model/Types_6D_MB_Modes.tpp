/*
 * Types_6DMB_Modes.tpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_6D_MB_MODES_TPP_
#define SRC_TYPES_6D_MB_MODES_TPP_

#include <ostream>
#include <iomanip>
#include <sstream>
#include "Types_6D_MB_Modes.h"

namespace as {

template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_6D_MB_Modes<REAL> const& dof)
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

	for (int lig = 0; lig < Common_Modes::numLigands; lig++){
		outStream 	<< setw(w) << dof._6D[lig].pos.x << setw(w) << dof._6D[lig].pos.y << setw(w) << dof._6D[lig].pos.z
					<< setw(w) << dof._6D[lig].ang.x << setw(w) << dof._6D[lig].ang.y << setw(w) << dof._6D[lig].ang.z;
		for(int mode=0;mode<Common_Modes::numModesLig[lig];mode++){
			outStream<< setw(w) << dof.modesLig[lig][mode];
		}
	}

	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result_6D_MB_Modes<REAL> const& enGrad) {
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


		s 	<< setw(width) << 0  << setw(width) << 0  << setw(width) << 0
		    << setw(width) << 0  << setw(width) << 0  << setw(width) << 0;
		for(int mode=0;mode<Common_Modes::numModesRec;mode++){
			s<< setw(width) << enGrad.modesRec[mode];
		}

		for (int lig = 0; lig < Common_Modes::numLigands; lig++){
			s << setw(width) << enGrad._6D[lig].ang.x  << setw(width) << enGrad._6D[lig].ang.y  << setw(width) << enGrad._6D[lig].ang.z
			  << setw(width) << enGrad._6D[lig].pos.x  << setw(width) << enGrad._6D[lig].pos.y  << setw(width) << enGrad._6D[lig].pos.z;

			for(int mode=0;mode<Common_Modes::numModesLig[lig];mode++){
				s<< setw(width) << enGrad.modesLig[lig][mode];
			}
		}
	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}


} // namespace as



#endif /* SRC_TYPES_6D_MB_TPP_ */
