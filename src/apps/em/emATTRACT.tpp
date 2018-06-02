/*
 * mcATTRACT.tpp
 *
 *  Created on: Jul 11, 2017
 *      Author: uwe
 */

#ifndef SRC_APPS_EM_EMATTRACT_TPP_
#define SRC_APPS_EM_EMATTRACT_TPP_

#include <iostream>
#include "emATTRACT.h"
#include "Configurator_6D.h"
#include "Request.h"
#include "Server.h"
#include "RequestHandler.h"



namespace as {

template<typename GenericTypes>
emATTRACT<GenericTypes>::emATTRACT() : _config(new configurator_t()) {}

template<typename GenericTypes>
void emATTRACT<GenericTypes>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename GenericTypes>
void emATTRACT<GenericTypes>::finalize() {
	_config->finalize();
}

template<typename GenericTypes>
void emATTRACT<GenericTypes>::run() {

	auto& dofs = _config->dofs();
	auto server = _config->server();
	auto& common = _config->common();

	RequestHandler<GenericTypes> requestHandler = RequestHandler<GenericTypes>::newBuilder()
			.withServer(server)
			.withCommon(common)
			.withDofs(dofs)
			.withSolverName("VA13")
			.build();

	requestHandler.run();

// THis is just for testing purposes
	auto results_grad = requestHandler.getResultEnGrads();
	auto results_dof = requestHandler.getResultStates();
	int count = 0;
	//for (auto const res : results_dof) {
	for( int i = 0; i< results_dof.size();i++){
		std::cout << results_grad[i].get_Energy() << std::endl;
		std::cout << results_dof[i] << std::endl;

	}

}

} //namespace



#endif /* SRC_APPS_EM_EMATTRACT_TPP_ */
