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
#include "Configurator_6D_Modes.h"
#include "Request.h"
#include "Server.h"
#include "RequestHandler.h"
//CHANGED
#include <fstream>
#include "Types_6D_Modes.h"
#include "Types_6D.h"
#include <string>
#include <sstream>
#include <ostream>


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


	bool use_modes = false ;

	if( std::is_same<GenericTypes, Types_6D_Modes<float>>::value  ||  std::is_same<GenericTypes, Types_6D_Modes<double>>::value){
		 use_modes = true;
	}
	std::string filename;


	auto dof_results = requestHandler.getResultEnGrads();
	auto results = requestHandler.getResultStates();


	for (auto const res : dof_results) {
		std::cout << res << std::endl;
	}

	std::string const proteinPath = "/home/glenn/Documents/Masterthesis/testfolder/2FD6";

	std::fstream fs;
	std::stringstream ss;
	if (use_modes){
		filename = proteinPath + "/output/GATTRACT/CPU/MODES/EM_OUT_PY.dat";
		ss << " " << "energy_total";
		ss <<" " << "rec_rot_x"<< " " << "rec_rot_y"<< " " << "rec_rot_z"<< " " << "rec_trans_x"<< " " << "rec_trans_y"<< " " << "rec_trans_z" ;
		ss << " " << "rec_mode_1"<< " " << "rec_mode_2"<< " " << "rec_mode_3"<< " " << "rec_mode_4"<< " " << "rec_mode_5"<< " ";
		ss << " " << "lig_rot_x"<< " " << "lig_rot_y"<< " " << "lig_rot_z"<< " " << "lig_trans_x"<< " " << "lig_trans_y"<< " " << "lig_trans_z";
		ss << " " << "lig_mode_1"<< " " << "lig_mode_2"<< " " << "lig_mode_3"<< " " << "lig_mode_4"<< " " << "lig_mode_5"<< " ";
		ss << std::endl;
	}
	else{
		filename = proteinPath + "/output/GATTRACT/CPU/NOMODES/EM_OUT_PY.dat";
		ss << " " << "energy_total";
		ss << " " << "lig_rot_x"<< " " << "lig_rot_y"<< " " << "lig_rot_z"<< " " << "lig_trans_x"<< " " << "lig_trans_y"<< " " << "lig_trans_z";
		ss << std::endl;

	}
	fs.open (filename, std::fstream::in | std::fstream::out | std::fstream::trunc );


	for ( int i= 0; i < results.size(); i++){
		print_results( ss, results[i], dof_results[i]);
	}
	fs << ss.str();
	fs.close();

}

} //namespace



#endif /* SRC_APPS_EM_EMATTRACT_TPP_ */
