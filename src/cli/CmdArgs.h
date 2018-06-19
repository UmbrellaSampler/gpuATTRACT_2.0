/*
 * CmdArgs.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_CMDARGS_H_
#define SRC_CMDARGS_H_

#include <vector>
#include <iosfwd>
#include <AppType.h>

namespace as {

class CmdArgs {
public:
	AppType app;

	/* Input Files */
	std::string dofName;
	std::string gridRecName;
	std::string gridLigName;
	std::string ligName;
	std::string recName;
	std::string paramsName;
	std::string alphabetRecName;
	std::string alphabetLigName;
	std::string recModesName;
	std::string ligModesName;
 
	/* Concurrency */
  int numModes;
	int numCPUs;
	std::vector<int> deviceIds;
	int chunkSize;

	/* Simulation */
	std::string precision;
	std::string dielec;
	double epsilon;
	double modeForceFac;
	double cutoff;

	/* Minimization */
	std::string solver;
	int maxConcurrency;
	int numChunks;

	friend std::ostream& operator<<(std::ostream& s, CmdArgs const& args);
};

}  // namespace as





#endif /* SRC_CMDARGS_H_ */
