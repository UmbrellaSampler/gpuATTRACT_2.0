/*
 * Configurator_6D_Modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_SERVICE_CONFIGURATOR_6D_MODES_H_
#define SRC_SERVICE_CONFIGURATOR_6D_MODES_H_

/*
 * Init_6D.h
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */


#include <memory>
#include "CmdArgs.h"
#include "Types_6D_modes.h"
#include "Service.h"

//template<typename GenericTypes>
//class Service;

namespace as {

template<typename GenericTypes>
class Server;

template<typename REAL>
class Configurator_6D_Modes {
	using service_t = Service<Types_6D_Modes<REAL>>;
//	using service_t = SERVICE;
	using input_t = typename service_t::input_t;
	using common_t = typename service_t::common_t;
//	using typename service_t::result_t;
public:
	using real_t = typename TypeWrapper<REAL>::real_t;
	using server_t = Server<Types_6D_Modes<REAL>>;

	std::vector<input_t>& dofs() noexcept {
		return _dofs;
	}

	std::shared_ptr<server_t> server() noexcept {
		return _server;
	}

	common_t& common() noexcept {
		return _ids;
	}

	void init(CmdArgs const& args) noexcept;
	void finalize() noexcept;
private:

	common_t _ids;
	std::vector<input_t> _dofs;
	std::shared_ptr<server_t> _server;

};

}  // namespace as




#endif /* CONFIGURATOR_6D_MODES_H_ */
