#ifndef SRC_REQUESTHANDLER_TPP_
#define SRC_REQUESTHANDLER_TPP_

#include <iostream>
#include <cassert>
#include <algorithm>

#include "nvToolsExt.h"

#include "RequestHandler.h"
#include "Chunk.h"
#include "SolverFactoryImpl.h"
#include "SolverBase.h"
#include "BFGSSolver.h"


//TODO: remove iostream
using std::cerr;
using std::cout;
using std::endl;

namespace as {

template<typename SERVER>
RequestHandler<SERVER>::RequestHandler(std::shared_ptr<extServer> server,
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
	init(solverName, dofs);
}

template<typename SERVER>
void RequestHandler<SERVER>::init(std::string const& solverName, std::vector<extDOF> const& dofs) {
	std::unique_ptr<SolverFactory> factory(new SolverFactoryImpl);

	/* Fill the object array */
	for (unsigned i = 0; i < dofs.size(); ++i) {
		SharedSolver solver;
		solver = factory->createSolverByName(solverName);

		solver->setState(dofs[i]);
		_objects.emplace_hint(_objects.end(), i, solver);
	}

	/* set number of recieved objects */
	_numObjects = _objects.size();

	/*
	 * initialize chunk list based on the number of available structures
	 */


	/* shrink the number of chunks to fit the minimal chunkSize */
	while (_numObjects < _numChunks*_minChunkSize && _numChunks > 1) {
		 --_numChunks;
	}
	assert(_numChunks >= 1);

	/* calculate chunk sizes */
	unsigned base_numPerChunk = MIN(_numObjects,_numConcurrentObjects) / _numChunks;
	unsigned rest = MIN(_numObjects,_numConcurrentObjects) % _numChunks;

	unsigned chunkSizes[_numChunks];
	std::fill(chunkSizes, chunkSizes + _numChunks, base_numPerChunk);

	assert(rest < _numChunks);
	for (unsigned i = 0; i < rest; ++i) {
		++chunkSizes[i];
	}

	/* setup chunks and put them into chunk list */
	for (unsigned i = 0; i < _numChunks; ++i) {
		_chunkList.emplace_back();
	}

	unsigned count = 0;
	for (auto& chunk : _chunkList) {
		ObjMapIter mapIter = _objects.begin();
		for (unsigned i = 0; i < chunkSizes[count]; ++i, ++mapIter) {
			chunk.getContainer().push_back(std::move(*mapIter));
			_objects.erase(mapIter);
			assert(mapIter != _objects.end());
		}
		assert(count < _numChunks);
		++count;

	}
}

//#define H_IO

template<typename SERVER>
void RequestHandler<SERVER>::run() {

	_collectedRequests.reserve(_chunkList.begin()->size());
	_collectedResults.reserve(_chunkList.begin()->size());

//	RingArray<int> reqIds;

	/* initial loop: start solvers and collect first requests and submit*/
	for (auto& chunk : _chunkList) {
		nvtxRangePushA("Processing");
		for (auto& obj : chunk.getContainer()) {
			SharedSolver& solver = obj.second;
			solver->start();
			_collectedRequests.push_back(TypesConverter<extDOF,Vector>::toFirst(solver->getState()));
		}
		nvtxRangePop();

		nvtxRangePushA("Submit");
//		int reqId = as::server_submit(*_server, _collectedRequests.data(), _collectedRequests.size(),
//				_serverOpt.gridId, _serverOpt.recId, _serverOpt.ligId, _serverOpt.useMode);

		request_t request(_collectedRequests.data(), _collectedRequests.size(), _common);
		try {
			_server->submit(request);
		} catch (std::invalid_argument &e) {
			cerr << "Error while submitting request: " << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}

		nvtxRangePop();
		chunk.setFetchSize(chunk.size());
		chunk.setRequest(request);

		_collectedRequests.resize(0);
	}


	unsigned count = 1;
	while(_finishedObjects.size () < _numObjects && count < 100000) {

		/* load balancing */

		/* Adjust chunk sizes to balance the workload of each chunk.
		 * This happens in case that the global object list is empty and
		 * the chunks cannot be refilled by new initial configurations.
		 * Do it not each iteration */
		if (_objects.empty() && count%4 == 0) {
			double ratio = Chunk::chunkSizeRatio(_chunkList);
			if(ratio > 1.5) {
				Chunk::loadBalanceChunks(_chunkList);
			}
		}


			for (auto chunkListIter = _chunkList.begin(); chunkListIter != _chunkList.end(); ) {
				auto& chunk = *chunkListIter;
				_collectedRequests.resize(0);
				if (chunk.fetchSize() > 0) {
					_collectedResults.resize(chunk.fetchSize());


					/* Wait for requests */
					nvtxRangePushA("Waiting");
					try {
						_server->wait(chunk.request(), _collectedResults.data());
					} catch (std::exception &e) {
						cerr << "Error while waiting for request to finish: " << e.what() << std::endl;
						exit(EXIT_FAILURE);
					}
					nvtxRangePop();

					/* Assigne results */
					chunk.setResults(_collectedResults);
				}

			/* Check if other chunks assigned results (after loadbalancing)*/
			chunk.checkLBconts();

			/* Process chunk and remove converged structures */
			nvtxRangePushA("Processing");
			auto iter = chunk.getContainer().begin();
			iter = chunk.getContainer().begin();

			while (iter != chunk.getContainer().end()) {
				SharedSolver& solver = iter->second;
				solver->step();

				/* test for convergence */
				if(solver->converged()) {

					/* destroy coroutine context by calling finalize */
					solver->finalize();
					/* move structure/object in finished object container */
					_finishedObjects.insert(move(*iter)); // no copy-construction
					chunk.getContainer().erase(iter++);

					/* move new structure/solver from object map if any left*/
					if (!_objects.empty()) {

						ObjMapIter objIter = _objects.begin();
						iter = chunk.getContainer().insert(iter, std::move(*objIter));
						_objects.erase(objIter);

						/* prepare new solver */
						SharedSolver& newSolver = iter->second;
						newSolver->start();
						/* collect new request */
						_collectedRequests.push_back(TypesConverter<extDOF,Vector>::toFirst(newSolver->getState()));
						++iter;
					}
				} else {
					/* collect new request */
					_collectedRequests.push_back(TypesConverter<extDOF,Vector>::toFirst(solver->getState()));
					++iter;
				}

			}
			nvtxRangePop();
			assert(iter == chunk.getContainer().end());

			chunk.setFetchSize(chunk.size());

			/* submit request */
			if (_collectedRequests.size() > 0) { // there is still something to submit
				nvtxRangePushA("Submit");
				request_t request(_collectedRequests.data(), _collectedRequests.size(), _common);
				try {
					_server->submit(request);
				} catch (std::invalid_argument &e) {
					cerr << "Error while submitting request: " << e.what() << std::endl;
					exit(EXIT_FAILURE);
				}

				nvtxRangePop();
				chunk.setRequest(request);
			}

			++chunkListIter;
			// do not remove any chunks since objects might still reside in LBconts waiting for results from other chunks

		} // for each chunk

		++count;
	} // while
	assert(_finishedObjects.size() == _numObjects);

}

template<typename SERVER>
auto RequestHandler<SERVER>::getResultStates() noexcept -> std::vector<extDOF>  {
	std::vector<extDOF> stateVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		stateVec[i] = TypesConverter<extDOF,Vector>::toFirst(_finishedObjects[i]->getState());
	}
	return stateVec;
}

template<typename SERVER>
auto RequestHandler<SERVER>::getResultEnGrads() noexcept -> std::vector<extEnGrad> {
	std::vector<extEnGrad> enGradVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		enGradVec[i] = TypesConverter<extEnGrad, ObjGrad>::toFirst(_finishedObjects[i]->getObjective());
	}
	return enGradVec;
}

template<typename SERVER>
auto RequestHandler<SERVER>::getStatistics() noexcept -> std::vector<std::unique_ptr<Statistic>> {
	std::vector<std::unique_ptr<Statistic>> statisticVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		statisticVec[i] = _finishedObjects[i]->getStats();
	}
	return statisticVec;
}

template<typename SERVER>
class RequestHandler<SERVER>::Builder {
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

} // namespace

#endif // SRC_REQUESTHANDLER_TPP_

