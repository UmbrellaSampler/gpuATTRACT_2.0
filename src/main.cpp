/*
 * main.cpp
 *
 *  Created on: Dec 29, 2015
 *      Author: uwe
 */


#include <iostream>
#include "CmdParser.h"
#include "CmdArgs.h"


using namespace std;
using namespace as;

int main(int argc, char* argv[]) {

	CmdParser parser;
	parser.parse(argc, argv);
	auto args = parser.args();
	std::cout << *args << std::endl;
	return 0;

}

