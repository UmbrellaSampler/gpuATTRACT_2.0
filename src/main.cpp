/*
 * main.cpp
 *
 *  Created on: Dec 29, 2015
 *      Author: uwe
 */


#include <iostream>
#include <exception>
#include <memory>
#include "CmdParser.h"
#include "CmdArgs.h"
#include "App.h"
#include "AppFactory.h"



using namespace as;

int main(int argc, char* argv[]) {

#ifdef CUDA
	std::cerr << "Using CUDA" << std::endl;
#else
	std::cerr << "Not Using CUDA" << std::endl;
#endif

	try {

		CmdParser parser;
		parser.parse(argc, argv);
		auto args = parser.args();
//		std::cout << args << std::endl;

		auto app = AppFactory::create(args);
		app->init(args);
		app->run();
	} catch (std::exception& e) {
		std::cerr << "Catched error in main: " << + e.what() << std::endl;
		exit(EXIT_FAILURE);
	}

	return 0;

}

