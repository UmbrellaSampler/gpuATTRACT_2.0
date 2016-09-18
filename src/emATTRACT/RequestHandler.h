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

#include <AttractServer>

namespace ema {

constexpr unsigned maxConcurrentObjects = 20000; // default maximum number of running coroutines that may exist at the same time.
constexpr unsigned numChunks = 2; // default number of chunks running at the same time. Each chunk maintains maxConcurrentObjects/numChunks objects.
constexpr unsigned minChunkSize = 10; // minimum chunksize that is worth to work with

template<typename SERVER>
class RequestHandler {

	using extDOF = typename SERVER::input_t;
	using extEnGrad = typename SERVER::result_t;
	using common_t = typename SERVER::common_t;
	using extServer = SERVER;
public:

	RequestHandler() : _server(nullptr),
		_numObjects(0), _numConcurrentObjects(maxConcurrentObjects), _numChunks(numChunks),
		_minChunkSize(minChunkSize){};

	void setNumConcurrentObjects(unsigned value) {_numConcurrentObjects = value;}
	void setNumChunks(unsigned value) {_numChunks = value;}
	void setMinChunkSize(unsigned value) { _minChunkSize = value;}

	/*
	 ** @brief: Initializes the RequestHandler. Member run() may now be called.
	 */
	void init(extServer& server, std::string const& solverName, std::vector<extDOF>& dofs, common_t common);

	void run();

	std::vector<extDOF> getResultStates();
	std::vector<extEnGrad> getResultEnGrads();
	std::vector<std::unique_ptr<Statistic>> getStatistics();

private:

	SERVER* _server;

	using SharedSolver = std::shared_ptr<SolverBase>;

	using ObjMap = std::map<unsigned, SharedSolver>;
	using ObjMapIter = ObjMap::iterator;
	ObjMap _objects;
	ObjMap _finishedObjects;

	using ChunkIter = Chunk::iterator;
	std::list<Chunk> _chunkList;

	std::vector<extDOF> _collectedRequests;
	std::vector<extEnGrad> _collectedResults;

	unsigned _numObjects;
	unsigned _numConcurrentObjects;
	unsigned _numChunks;
	unsigned _minChunkSize;
};

} // namespace


#endif /* REQUESTHANDLER_H_ */
