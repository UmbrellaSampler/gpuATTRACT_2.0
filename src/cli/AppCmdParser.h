/*
 * AppCmdParser.h
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#ifndef SRC_APPCMDPARSER_H_
#define SRC_APPCMDPARSER_H_

#include <memory>
#include <boost/program_options.hpp>
#include "CmdArgs.h"

namespace as {

class AppCmdParser {
public:
	AppCmdParser(std::shared_ptr<CmdArgs> args) : _args(args) {}
	virtual ~AppCmdParser() {};

	void parse(int argc, char* argv[]);

protected:
	virtual std::string appShortName() const noexcept = 0;

	virtual void addOptions() noexcept = 0;
	virtual void enforceRules(boost::program_options::variables_map const& vm) const = 0;
	virtual void assignArgs(boost::program_options::variables_map const& vm) noexcept = 0;

private:
	std::string usage() const noexcept;

	void printHelp() const noexcept;

protected:
	std::shared_ptr<CmdArgs> _args;
	boost::program_options::options_description _optsDesc;

private:
	std::string _applicationName;

};

} // namespace



#endif /* SRC_APPCMDPARSER_H_ */
