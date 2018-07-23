/*
 * mcATTRACT.tpp
 *
 *  Created on: Jul 11, 2017
 *      Author: uwe
 */

#ifndef SRC_APPS_MC_MCATTRACT_TPP_
#define SRC_APPS_MC_MCATTRACT_TPP_

#include <iostream>
#include "mcATTRACT.h"
#include "Configurator_6D.h"
#include "Request.h"
#include "Server.h"
#include "RequestHandler.h"
#include <chrono>



namespace as {

template<typename GenericTypes>
mcATTRACT<GenericTypes>::mcATTRACT() : _config(new configurator_t()) {}

template<typename GenericTypes>
void mcATTRACT<GenericTypes>::init(CmdArgs const& args) {
	_config->init(args);
}

template<typename GenericTypes>
void mcATTRACT<GenericTypes>::finalize() {
	_config->finalize();
}

template<typename GenericTypes>
void mcATTRACT<GenericTypes>::run() {

	auto& dofs = _config->dofs();
	auto server = _config->server();
	auto& common = _config->common();

		RequestHandler<GenericTypes> requestHandler = RequestHandler<GenericTypes>::newBuilder()
				.withServer(server)
				.withCommon(common)
				.withDofs(dofs)
				.withSolverName("MC")
				.build();
	auto start = std::chrono::system_clock::now();
	requestHandler.run();
	auto end = std::chrono::system_clock::now();
// THis is just for testing purposes
	auto results_grad = requestHandler.getResultEnGrads();
	auto results_dof = requestHandler.getResultStates();
	int count = 0;
	//for (auto const res : results_dof) {


//	std::cout << "#pivot 1 "<< " "<< 8.95646 << " "<<56.9144 << " "<<92.5237 << std::endl;
//	std::cout << "#pivot 2 "<< " "<<8.91543 << " "<<34.5702 << " "<<114.345 << std::endl;

	std::cout << "#pivot 1 "<< " "<< common.pivotRec.x << " "<<common.pivotRec.y << " "<<common.pivotRec.z << std::endl;
	std::cout << "#pivot 2 "<< " "<<common.pivotLig.x << " "<<common.pivotLig.y << " "<<common.pivotLig.z << std::endl;
	std::cout << "#centered receptor: false "<< std::endl;
	std::cout << "#centered ligand: false "<< std::endl;
	for (int i = 0; i < results_dof.size(); i++){
		std::cout << "#"<< i+1<< std::endl;
		std::cout <<"## Energy: " << results_grad[i].get_Energy() << std::endl;
		auto  pos = results_dof[i].get_pos();
				pos.x = pos.x + common.pivotRec.x -  common.pivotLig.x;
				pos.y = pos.y + common.pivotRec.y -  common.pivotLig.y;
				pos.z = pos.z + common.pivotRec.z -  common.pivotLig.z;
				results_dof[i].set_pos(pos.x,pos.y, pos.z);
				std::cout << results_dof[i] << std::endl;


	}


		//std::cout << "time"<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
}

} //namespace



#endif /* SRC_APPS_MC_MCATTRACT_TPP_ */
