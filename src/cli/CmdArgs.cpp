/*
 * CmdArgs.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <ostream>

#include "CmdArgs.h"

using namespace std;

namespace as {

std::ostream& operator<<(std::ostream& s, CmdArgs const& args) {
	s << "dofName     = "  << args.dofName             << endl;
	s << "gridName    = "  << args.gridName            << endl;
	s << "ligName     = "  << args.ligName             << endl;
	s << "recName     = "  << args.recName             << endl;
	s << "paramsName  = "  << args.paramsName          << endl;
	s << "alphabeName = "  << args.alphabetName 	   << endl;
	s << "numCPUs     = "  << args.numCPUs             << endl;
	s << "devices     = [ "; for (auto device : args.deviceIds) s << device << " "; s << "]"<< endl;
	s << "chunkSize   = "  << args.chunkSize           << endl;
	s << "prec        = "  << args.precision		   << endl;
	s << "dielec      = "  << args.dielec              << endl;
	s << "epsilon     = "  << args.epsilon;
	return s;
}

}  // namespace as


