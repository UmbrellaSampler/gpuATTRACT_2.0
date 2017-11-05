/*
 * AppParserFactory.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <exception>
#include "AppParserFactory.h"
#include "ScoreParser.h"
#include "EmParser.h"

using namespace as;

std::unique_ptr<AppCmdParser> AppParserFactory::create(AppType app, std::shared_ptr<CmdArgs> args) {
	std::unique_ptr<AppCmdParser> parser;

	switch(app) {
	case AppType::SCORE:
		parser = std::unique_ptr<ScoreParser>(new ScoreParser(args));
		break;
	case AppType::EM:
		parser = std::unique_ptr<EmParser>(new EmParser(args));
		break;
	default:
		throw std::invalid_argument("unknown app to create: " + static_cast<int>(app));
	}

	return parser;
}




