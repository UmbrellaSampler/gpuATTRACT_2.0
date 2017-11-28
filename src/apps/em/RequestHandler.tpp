#ifndef SRC_REQUESTHANDLER_TPP_
#define SRC_REQUESTHANDLER_TPP_

#include <iostream>
#include <cassert>
#include <algorithm>

#include "Chunk.h"
#include "RequestHandler.h"

#include "SolverFactoryImpl.h"
#include "SolverBase.h"
#include "BFGSSolver.h"


//TODO: remove iostream
using std::cerr;
using std::cout;
using std::endl;

namespace as {

template<typename GenericTypes>
RequestHandler<GenericTypes>::RequestHandler(std::shared_ptr<extServer> server,
			unsigned numConcurrentObjects,
			unsigned numChunks,
			unsigned minChunkSize,
			std::vector<extDOF> const& dofsInput,
			common_t const& common,
			std::string const& solverName) :
		_server(server),
		_numConcurrentObjects(numConcurrentObjects),
		_numChunks(numChunks),
		_minChunkSize(minChunkSize),
		_numObjects(0),
		_common(common)
{

	init(solverName, dofsInput);
}

template<typename GenericTypes>
void RequestHandler<GenericTypes>::init(std::string const& solverName, std::vector<extDOF> const& dofs) {
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

template<typename GenericTypes>
void RequestHandler<GenericTypes>::run() {

	_collectedRequests.reserve(_chunkList.begin()->size());
	_collectedResults.reserve(_chunkList.begin()->size());

//	RingArray<int> reqIds;

	/* initial loop: start solvers and collect first requests and submit*/
	for (auto& chunk : _chunkList) {
		for (auto& obj : chunk.getContainer()) {
			SharedSolver& solver = obj.second;
			solver->start();
			extDOF dof;
			dof.setDOFfromVector(solver->getState(),solver->getInputDOFConfig());
			_collectedRequests.push_back(dof);
			//_collectedRequests.push_back(TypesConverter<extDOF,Vector>::toFirst(solver->getState()));
		}

		request_t request(_collectedRequests.data(), _collectedRequests.size(), _common);
		try {
			_server->submit(request);
		} catch (std::invalid_argument &e) {
			cerr << "Error while submitting request: " << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}

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
					try {
						_server->wait(chunk.request(), _collectedResults.data());
					} catch (std::exception &e) {
						cerr << "Error while waiting for request to finish: " << e.what() << std::endl;
						exit(EXIT_FAILURE);
					}

					/* Assigne results */
					chunk.setResults(_collectedResults);
				}

			/* Check if other chunks assigned results (after loadbalancing)*/
			chunk.checkLBconts();

			/* Process chunk and remove converged structures */
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
						extDOF dof;
						dof.setDOFfromVector(newSolver->getState(),newSolver->getInputDOFConfig());
						_collectedRequests.push_back(dof);
						//_collectedRequests.push_back(TypesConverter<extDOF,Vector>::toFirst(newSolver->getState()));
						++iter;
					}
				} else {
					/* collect new request */
					extDOF dof;
					dof.setDOFfromVector(solver->getState(),solver->getInputDOFConfig());
					_collectedRequests.push_back(dof);
				//	_collectedRequests.push_back(TypesConverter<extDOF,Vector>::toFirst(solver->getState()));
					++iter;
				}

			}
			assert(iter == chunk.getContainer().end());

			chunk.setFetchSize(chunk.size());

			/* submit request */
			if (_collectedRequests.size() > 0) { // there is still something to submit
				request_t request(_collectedRequests.data(), _collectedRequests.size(), _common);
				try {
					_server->submit(request);
				} catch (std::invalid_argument &e) {
					cerr << "Error while submitting request: " << e.what() << std::endl;
					exit(EXIT_FAILURE);
				}
				chunk.setRequest(request);
			}

			++chunkListIter;
			// do not remove any chunks since objects might still reside in LBconts waiting for results from other chunks

		} // for each chunk

		++count;
	} // while
	assert(_finishedObjects.size() == _numObjects);

}

template<typename GenericTypes>
auto RequestHandler<GenericTypes>::getResultStates() noexcept -> std::vector<extDOF>  {
	std::vector<extDOF> stateVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		//stateVec[i] = TypesConverter<extDOF,Vector>::toFirst(_finishedObjects[i]->getState());
		extDOF dof;
		dof.setDOFfromVector(_finishedObjects[i]->getState(),_finishedObjects[i]->getInputDOFConfig());
		stateVec[i] = dof;
	}
	return stateVec;
}

template<typename GenericTypes>
auto RequestHandler<GenericTypes>::getResultEnGrads() noexcept -> std::vector<extEnGrad> {
	std::vector<extEnGrad> enGradVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		extEnGrad dof;
		dof.setDOFfromVector(_finishedObjects[i]->getObjective(),_finishedObjects[i]->getResultDOFConfig());
		enGradVec[i] = dof;
		//enGradVec[i] = TypesConverter<extEnGrad, ObjGrad>::toFirst(_finishedObjects[i]->getObjective());
	}
	return enGradVec;
}

template<typename GenericTypes>
auto RequestHandler<GenericTypes>::getStatistics() noexcept -> std::vector<std::unique_ptr<Statistic>> {
	std::vector<std::unique_ptr<Statistic>> statisticVec(_finishedObjects.size());
	for (unsigned i = 0; i < _finishedObjects.size(); ++i) {
		statisticVec[i] = _finishedObjects[i]->getStats();
	}
	return statisticVec;
}


} // namespace

#endif // SRC_REQUESTHANDLER_TPP_

