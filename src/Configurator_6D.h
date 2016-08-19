/*
 * Init_6D.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#ifndef SRC_INIT_6D_H_
#define SRC_INIT_6D_H_

#include <memory>
#include "CmdArgs.h"
#include "CPU_6D_EnergyService.h"

namespace as {

template<typename SERVICE>
class Server;

template<typename REAL>
class Configurator_6D {
	using service_t = CPU_6D_EnergyService<REAL>;
	using dof_t = typename service_t::dof_t;
	using common_t = typename service_t::common_t;
	using result_t = typename service_t::result_t;
public:
	using server_t = Server<service_t>;

	std::vector<dof_t>& dofs() {
		return _dofs;
	}

	server_t& server() {
		return *_server;
	}

	common_t& common() {
		return _ids;
	}

	void init(CmdArgs const& args);
	void finalize();
private:

	common_t _ids;
	std::vector<dof_t> _dofs;
	std::shared_ptr<server_t> _server;

};

}  // namespace as



#endif /* SRC_INIT_6D_H_ */
