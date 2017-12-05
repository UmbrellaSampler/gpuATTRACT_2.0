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
	Request<input_t, common_t> request(dofs.data(), numDofs, common);
	server->submit(request);

	auto results = std::vector<result_t>(dofs.size());
	server->wait(request, results.data());

	for (result_t const res : results) {
		std::cout << res << std::endl;
	}

}

}  // namespace as



#endif /* SRC_SCATTRACT_TPP_ */
