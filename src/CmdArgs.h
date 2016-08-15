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

namespace as {

enum class AppType {
	Score
};

class CmdArgs {
public:
	AppType app;
	std::string dofName;
	std::string gridName;
	std::string ligName;
	std::string recName;
	std::string paramsName;
	std::string recGridAlphabetName;
	int numCPUs;
	std::vector<int> deviceIds;
	int chunkSize;
	std::string precision;

	friend std::ostream& operator<<(std::ostream& s, CmdArgs const& args);
};

}  // namespace as





#endif /* SRC_CMDARGS_H_ */
