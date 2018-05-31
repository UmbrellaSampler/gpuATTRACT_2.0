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


	auto results = requestHandler.getResultEnGrads();

	for (result_t const res : results) {
		std::cout << res << std::endl;
	}

}

} //namespace



#endif /* SRC_APPS_EM_EMATTRACT_TPP_ */
