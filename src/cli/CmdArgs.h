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
#include "Types_6D_Config.h"
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

	std::string alphabetName;
	std::string modeNames[NUM_MAX_PROTEIN];
	std::string proteinNames[NUM_MAX_PROTEIN];
	std::string gridNames[NUM_MAX_PROTEIN];
 
	/* Concurrency */
  int numModes;
	int numCPUs;
	std::vector<int> deviceIds;
	int chunkSize;
	int numProtein;

	/* Simulation */
	std::string precision;
	std::string dielec;
	double epsilon;
	double modeEVFactor;
	double radius_cutoff;

	/* Minimization */
	std::string solver;
	int maxConcurrency;
	int numChunks;

	friend std::ostream& operator<<(std::ostream& s, CmdArgs const& args);
};

}  // namespace as





#endif /* SRC_CMDARGS_H_ */
