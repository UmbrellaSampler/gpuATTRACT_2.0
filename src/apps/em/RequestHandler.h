/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef REQUESTHANDLER_H_
#define REQUESTHANDLER_H_

#include <vector>
#include <map>
#include <list>
#include <memory>

namespace as {

class SolverBase;
struct Statistic;

template<typename SERVER>
class RequestHandler {
	class Builder;
	class Chunk;

public:
	static constexpr unsigned DEFAULT_MAX_CONCURRENT_OBJECTS = 20000; // default maximum number of asynchronous coroutines
	static constexpr unsigned DEFAULT_NUM_CHUNKS = 2; // default number of chunks running at the same time. Each chunk maintains numConcurrentObjects/numChunks objects.
	static constexpr unsigned DEFAULT_MIN_CHUNK_SIZE = 10; // minimum chunksize that is worth to work with

private:
	using extDOF = typename SERVER::input_t;
	using extEnGrad = typename SERVER::result_t;
	using common_t = typename SERVER::common_t;
	using extServer = SERVER;
	using SharedSolver = std::shared_ptr<SolverBase>;
	using ObjMap = std::map<unsigned, SharedSolver>;
	using ObjMapIter = ObjMap::iterator;
//	using ChunkIter = typename Chunk::iterator;

	using request_t = typename SERVER::request_t;

	std::shared_ptr<extServer> _server;
	const unsigned _numConcurrentObjects;
	unsigned _numChunks;
	const unsigned _minChunkSize;

	unsigned _numObjects;
	ObjMap _objects;
	ObjMap _finishedObjects;

	const common_t _common;
	std::list<Chunk> _chunkList;
	std::vector<extDOF> _collectedRequests;
	std::vector<extEnGrad> _collectedResults;


	explicit RequestHandler(std::shared_ptr<extServer> server,
			unsigned numConcurrentObjects,
			unsigned numChunks,
			unsigned minChunkSize,
			std::vector<extDOF> const& dofs,
			common_t const& common,
			std::string const& solverName);

	void init(std::string const& solverName, std::vector<extDOF> const& dofs);

public:

	void run();

	std::vector<extDOF> getResultStates() noexcept;
	std::vector<extEnGrad> getResultEnGrads() noexcept;
	std::vector<std::unique_ptr<Statistic>> getStatistics() noexcept;

	static Builder newBuilder() {
		return Builder();
	}
};

} // namespace

#include "RequestHandlerBuilder.h"


#endif /* REQUESTHANDLER_H_ */
