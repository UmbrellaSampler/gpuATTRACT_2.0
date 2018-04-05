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
#include <memory>
#include <fstream>
#include <sstream>

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
	std::shared_ptr<Request<input_t, common_t>> sh_request = std::make_shared<Request<input_t, common_t>>(request);
	server->submit(sh_request);

	auto results = std::vector<result_t>(dofs.size());
	server->wait(sh_request, results.data());


	for (result_t const res : results) {
		std::cout << res << std::endl;
	}

	bool use_modes = false ;

	if( std::is_same<GenericTypes, Types_6D_Modes<float>>::value  ||  std::is_same<GenericTypes, Types_6D_Modes<double>>::value){
		 use_modes = true;
	}
	std::string filename;



	for (auto const res : results) {
		std::cout << res << std::endl;
	}

	std::fstream fs;
	std::stringstream ss;

	if (use_modes){
		filename = "/home/glenn/Documents/Masterthesis/testfolder/1AVX/output/GATTRACT/CPU/MODES/SCORE_OUT_PY.dat";
		ss << " " << "energy_total";
		ss <<" " << "rec_rot_x"<< " " << "rec_rot_y"<< " " << "rec_rot_z"<< " " << "rec_trans_x"<< " " << "rec_trans_y"<< " " << "rec_trans_z" ;
		ss << " " << "rec_mode_1"<< " " << "rec_mode_2"<< " " << "rec_mode_3"<< " " << "rec_mode_4"<< " " << "rec_mode_5"<< " ";
		ss << " " << "lig_rot_x"<< " " << "lig_rot_y"<< " " << "lig_rot_z"<< " " << "lig_trans_x"<< " " << "lig_trans_y"<< " " << "lig_trans_z";
		ss << " " << "lig_mode_1"<< " " << "lig_mode_2"<< " " << "lig_mode_3"<< " " << "lig_mode_4"<< " " << "lig_mode_5"<< " ";
		ss << std::endl;
	}
	else{
		filename = "/home/glenn/Documents/Masterthesis/testfolder/1AVX/output/GATTRACT/CPU/NOMODES/SCORE_OUT_PY.dat";
		ss << " " << "energy_total";
		ss << " " << "lig_rot_x"<< " " << "lig_rot_y"<< " " << "lig_rot_z"<< " " << "lig_trans_x"<< " " << "lig_trans_y"<< " " << "lig_trans_z";
		ss << std::endl;
	}
	fs.open (filename, std::fstream::in | std::fstream::out | std::fstream::trunc );

	for ( int i= 0; i < results.size(); i++){
		print_results( ss, results[i]);
	}

	fs << ss.str();
	fs.close();

}

}  // namespace as



#endif /* SRC_SCATTRACT_TPP_ */
