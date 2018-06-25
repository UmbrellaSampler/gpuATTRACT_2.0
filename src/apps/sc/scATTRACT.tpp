/*
 * scATTRACT.tpp
 *
 *  Created on: Aug 17, 2016
 *      Author: uwe
 */

#ifndef SRC_SCATTRACT_TPP_
#define SRC_SCATTRACT_TPP_

#include <iostream>
#include "scATTRACT.h"
#include "Configurator_6D.h"
#include "Request.h"
#include "Server.h"
#include <chrono>
#define DEBUG
#ifdef DEBUG
#include <fstream>
#include <sstream>
#include "debug_functions.h"
#include <string>
#include "ServiceFactory.h"
#endif
namespace as {

template<typename GenericTypes>
scATTRACT<GenericTypes>::scATTRACT() : _config(new configurator_t()) {}

template<typename GenericTypes>
void scATTRACT<GenericTypes>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename GenericTypes>
void scATTRACT<GenericTypes>::finalize() {
	_config->finalize();
}

template<typename GenericTypes>
void scATTRACT<GenericTypes>::run() {

	auto& dofs = _config->dofs();
	auto server = _config->server();
	auto& common = _config->common();
	size_t numDofs = dofs.size();
	auto request = std::make_shared<Request<input_t, common_t>>( Request<input_t, common_t>(dofs.data(), numDofs, common) );
	auto start = std::chrono::system_clock::now();

	/* do some work */


	server->submit(request);

	auto results = std::vector<result_t>(dofs.size());
	server->wait(request, results.data());
	auto end = std::chrono::system_clock::now();
	#ifdef DEBUG //_config->server()::service_t CPUEnergyService6DModes<real_t>
		bool  use_GPU = false;
			#ifdef CUDA
			use_GPU = true;
			#endif

		std::string file_debug = getDebugPath<real_t, result_t>( use_GPU );
		std::fstream fs;
		fs.open (file_debug, std::fstream::in | std::fstream::out | std::fstream::trunc );
		std::stringstream ss;
		print_header<real_t, result_t>( ss);
		ss << std::endl;
		for (result_t const res : results)
		{
			print_enGrad<real_t, result_t>( ss, res);
			ss << std::endl;


		}
		fs << ss.str();
		fs.close();
	#endif

	std::cout << "time"<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
	for (result_t const res : results) {
		std::cout << res << std::endl;
	}

}

}  // namespace as



#endif /* SRC_SCATTRACT_TPP_ */
