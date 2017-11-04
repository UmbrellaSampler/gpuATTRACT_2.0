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
	std::string dofName;
	std::string gridName;
	std::string ligName;
	std::string recName;
	std::string paramsName;
	std::string alphabetName;
	int numModes;
	std::string recModesName;
	std::string ligModesName;
	int numCPUs;
	std::vector<int> deviceIds;
	int chunkSize;
	std::string precision;
	std::string dielec;
	double epsilon;

	friend std::ostream& operator<<(std::ostream& s, CmdArgs const& args);
};

}  // namespace as





#endif /* SRC_CMDARGS_H_ */
