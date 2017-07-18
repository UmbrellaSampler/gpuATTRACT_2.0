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

namespace as {

template<typename SERVICE>
class Server;

template<typename SERVICE>
class Configurator_6D {
	using service_t = SERVICE;
	using real_t = typename service_t::real_t;
	using dof_t = typename service_t::dof_t;
	using common_t = typename service_t::common_t;
public:
	using server_t = Server<service_t>;

	std::vector<dof_t>& dofs() {
		return _dofs;
	}

	std::shared_ptr<server_t> server() {
		return _server;
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
