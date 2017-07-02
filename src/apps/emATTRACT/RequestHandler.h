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
#include "SolverBase.h"
#include "BFGSSolver.h"
#include "Chunk.h"

namespace as {


template<typename SERVER>
class RequestHandler {
public:
	constexpr unsigned DEFAULT_MAX_CONCURRENT_OBJECTS = 20000; // default maximum number of asynchronous coroutines
	constexpr unsigned DEFAULT_NUM_CHUNKS = 2; // default number of chunks running at the same time. Each chunk maintains numConcurrentObjects/numChunks objects.
	constexpr unsigned DEFAULT_MIN_CHUNK_SIZE = 10; // minimum chunksize that is worth to work with

private:
	using extDOF = typename SERVER::input_t;
	using extEnGrad = typename SERVER::result_t;
	using common_t = typename SERVER::common_t;
	using extServer = SERVER;
	using SharedSolver = std::shared_ptr<SolverBase>;
	using ObjMap = std::map<unsigned, SharedSolver>;
	using ObjMapIter = ObjMap::iterator;
	using ChunkIter = Chunk::iterator;

	std::shared_ptr<extServer> _server;
	const unsigned _numConcurrentObjects;
	const unsigned _numChunks;
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
			std::string const& solverName) :
		_server(server),
		_numConcurrentObjects(numConcurrentObjects),
		_numChunks(numChunks),
		_minChunkSize(minChunkSize),
		_numObjects(0),
		_common(common)
	{
		init(solverName, dofs, common);
	};

public:

	void init(std::string const& solverName, std::vector<extDOF>& dofs, common_t const& common);

	void run();

	std::vector<extDOF> getResultStates();
	std::vector<extEnGrad> getResultEnGrads();
	std::vector<std::unique_ptr<Statistic>> getStatistics();



	class Builder {
	private:
		std::shared_ptr<extServer> _server;
		unsigned _numConcurrentObjects = DEFAULT_MAX_CONCURRENT_OBJECTS;
		unsigned _numChunks = DEFAULT_NUM_CHUNKS;
		unsigned _minChunkSize = DEFAULT_MIN_CHUNK_SIZE;

	public:
		Builder& withServer(std::shared_ptr<extServer> server) {
			_server = server;
			return *this;
		}

		Builder& withNumConcurrentObjects(unsigned numConcurrentObjects) {
			_numConcurrentObjects = numConcurrentObjects;
			return *this;
		}

		Builder& withNumChunks(unsigned numChunks) {
			_numChunks = numChunks;
			return *this;
		}

		Builder& withMinChunkSize(unsigned minChunkSize) {
			_minChunkSize = minChunkSize;
			return *this;
		}
	};

	static Builder newBuilder() {
		return Builder();
	}


};

} // namespace


#endif /* REQUESTHANDLER_H_ */
