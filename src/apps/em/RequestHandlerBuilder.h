/*
 * RequestHandlerBuilder.h
 *
 *  Created on: Jul 18, 2017
 *      Author: uwe
 */

#ifndef SRC_APPS_EM_REQUESTHANDLERBUILDER_H_
#define SRC_APPS_EM_REQUESTHANDLERBUILDER_H_

#include "RequestHandler.h"

namespace as {

template<typename SERVER>
class RequestHandler<SERVER>::Builder {
private:
	std::shared_ptr<extServer> _server;
	unsigned _numConcurrentObjects = DEFAULT_MAX_CONCURRENT_OBJECTS;
	unsigned _numChunks = DEFAULT_NUM_CHUNKS;
	unsigned _minChunkSize = DEFAULT_MIN_CHUNK_SIZE;
	std::vector<extDOF> const* _dofs;
	common_t _common;
	std::string _solverName;

public:
	Builder& withServer(std::shared_ptr<extServer> server) noexcept {
		_server = server;
		return *this;
	}

	Builder& withNumConcurrentObjects(unsigned numConcurrentObjects) noexcept {
		_numConcurrentObjects = numConcurrentObjects;
		return *this;
	}

	Builder& withNumChunks(unsigned numChunks) noexcept {
		_numChunks = numChunks;
		return *this;
	}

	Builder& withMinChunkSize(unsigned minChunkSize) noexcept {
		_minChunkSize = minChunkSize;
		return *this;
	}

	Builder& withDofs(std::vector<extDOF> const& dofs) noexcept {
		_dofs = &dofs;
		return *this;
	}

	Builder& withCommon(common_t const& common) noexcept {
		_common = common;
		return *this;
	}

	Builder& withSolverName(std::string const& solverName) noexcept {
		_solverName = solverName;
		return *this;
	}


	RequestHandler build() noexcept {
		return RequestHandler(_server, _numConcurrentObjects,
				_numChunks, _minChunkSize,
				*_dofs, _common, _solverName);
	}

};

} // namespace

#endif /* SRC_APPS_EM_REQUESTHANDLERBUILDER_H_ */
