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

template<typename SERVICE>
emATTRACT<SERVICE>::emATTRACT() : _config(new configurator_t()) {}

template<typename SERVICE>
void emATTRACT<SERVICE>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename SERVICE>
void emATTRACT<SERVICE>::finalize() {
	_config->finalize();
}

template<typename SERVICE>
void emATTRACT<SERVICE>::run() {

	auto& dofs = _config->dofs();
	auto server = _config->server();
	auto& common = _config->common();

	RequestHandler<server_t> requestHandler = RequestHandler<server_t>::newBuilder()
			.withServer(server)
			.withCommon(common)
			.withDofs(dofs)
			.withSolverName("sdf")
			.build();

	requestHandler.run();


	auto results = requestHandler.getResultEnGrads();

	for (result_t const res : results) {
		std::cout << res << std::endl;
	}

}

} //namespace



#endif /* SRC_APPS_EM_EMATTRACT_TPP_ */
