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
#include "Configurator_6D_Modes.h"
#include "Request.h"
#include "Server.h"

namespace as {

template<typename SERVICE>
scATTRACT<SERVICE>::scATTRACT() : _config(new configurator_t()) {}

template<typename SERVICE>
void scATTRACT<SERVICE>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename SERVICE>
void scATTRACT<SERVICE>::finalize() {
	_config->finalize();
}

template<typename SERVICE>
void scATTRACT<SERVICE>::run() {

	auto& dofs = _config->dofs();
	auto server = _config->server();
	auto& common = _config->common();
	size_t numDofs = dofs.size();
	Request<input_t, common_t> request(dofs.data(), numDofs, common);
	server->submit(request);

	auto results = std::vector<result_t>(dofs.size());
	server->wait(request, results.data());

	for (result_t const res : results) {
	//	std::cout << res << std::endl;
	}

}

}  // namespace as



#endif /* SRC_SCATTRACT_TPP_ */
