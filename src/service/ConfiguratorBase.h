/*
 * ConfiguratorBase.h
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_SERVICE_CONFIGURATORBASE_H_
#define SRC_SERVICE_CONFIGURATORBASE_H_

#include <memory>
#include "CmdArgs.h"
#include "Types_6D.h"
#include "Service.h"

namespace as {

template<typename GenericTypes>
class Server;

template<typename GenericTypes>
class ConfiguratorBase {

protected:
	using service_t = Service<GenericTypes>;
	using input_t = typename service_t::input_t;
	using common_t = typename service_t::common_t;

public:
	using real_t = typename input_t::real_t;
	using server_t = Server<GenericTypes>;

	std::vector<input_t>& dofs() noexcept {
		return _dofs;
	}

	std::shared_ptr<server_t> server() noexcept {
		return _server;
	}

	common_t& common() noexcept {
		return _ids;
	}

	virtual void init(CmdArgs const& args) = 0;
	virtual void finalize() noexcept = 0;

	virtual ~ConfiguratorBase() {}

protected:
	common_t _ids;
	std::vector<input_t> _dofs;
	std::shared_ptr<server_t> _server;

};

}  // namespace as



#endif /* SRC_SERVICE_CONFIGURATORBASE_H_ */
