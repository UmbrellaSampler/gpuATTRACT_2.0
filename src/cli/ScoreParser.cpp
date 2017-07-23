/*
 * ScoreParser.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: uwe
 */

#include <sstream>

#include "ScoreParser.h"
#include "parser_constants.h"


using namespace as;
using namespace std;
namespace po = boost::program_options;

namespace as {

std::string ScoreParser::appShortName() const noexcept {
	return APP_SHORT_NAME_SC;
}

}
