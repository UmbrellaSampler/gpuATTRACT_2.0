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
#include "time.h"


using namespace as;

int main(int argc, char* argv[]) {
	clock_t start = clock();
#ifdef CUDA
	//std::cerr << "CUDA enabled" << std::endl;
#else
	std::cerr << "CUDA not enabled" << std::endl;
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
		std::cerr << "Error in main: " << + e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
	clock_t end = clock();
		double elapsed_time = (end - start)/(double)CLOCKS_PER_SEC;
		printf("elapsed total time: %f",elapsed_time);
	return 0;

}

