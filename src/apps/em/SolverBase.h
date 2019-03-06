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

#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_

#include <boost/coroutine/all.hpp>
#include <boost/version.hpp>
#include <Eigen/Core>
#include <cassert>
#include<Eigen/StdVector>
#include <meta.h>

namespace as {


#if BOOST_VERSION <= 105400
using coro_t = boost::coroutines::coroutine<void(void)>;
using pull_type = coro_t;
using push_type = coro_t::caller_type;
#elif BOOST_VERSION == 105500
using coro_t = boost::coroutines::coroutine<void>;
using pull_type = coro_t::pull_type;
using push_type = coro_t::push_type;
#elif BOOST_VERSION >= 105600
using coro_t = boost::coroutines::asymmetric_coroutine<void>;
using pull_type = coro_t::pull_type;
using push_type = coro_t::push_type;
#endif

//static_assert(false, BOOST_LIB_VERSION );

struct Statistic {
	virtual Statistic* getCopy() const = 0;
	friend std::ostream& operator << (std::ostream& os, const Statistic& stats) {
      return stats.print(os); // polymorphic print via reference
    }
	virtual ~Statistic() {};

	unsigned numRequests = 0;
private:
	virtual std::ostream& print(std::ostream&) const = 0;

};

class SolverBase {
public:
	SolverBase() : coro(nullptr),trackedStates(std::make_shared<std::vector<std::vector<float>>>()) ,trackedGrads(std::make_shared<std::vector<std::vector<float>>>()){}
	virtual ~SolverBase() { delete coro;}

	/* make object not copyable, but movealble only */
	SolverBase(const SolverBase& ) = delete;
	SolverBase& operator= (const SolverBase& ) = delete;

	SolverBase(SolverBase && rhs) {
		state = std::move(rhs.state);
		objective = std::move(rhs.objective);
		coro = std::move(rhs.coro);
		rhs.coro = nullptr;
		trackedStates= std::make_shared<std::vector<std::vector<float>>>();
		trackedGrads= std::make_shared<std::vector<std::vector<float>>>();
	}

	SolverBase& operator= (SolverBase&& rhs) {
		state = std::move(rhs.state);
		objective = std::move(rhs.objective);
		coro = std::move(rhs.coro);
		rhs.coro = nullptr;
		trackedStates= std::make_shared<std::vector<std::vector<float>>>();
		trackedGrads= std::make_shared<std::vector<std::vector<float>>>();

		return *this;
	}

	bool converged() {return !*coro;}
	void setState(const Vector& value) { state = value;}

	template <typename DOFType>
	void setState(const DOFType& value) { state = TypesConverter<DOFType, Vector>::toSecond(value);}

	Vector getState() {return state;}

	void setObjective(const ObjGrad& value) { objective = value; }

	template <typename ResultType>
	void setObjective(const ResultType& value) {	objective = TypesConverter<ResultType, ObjGrad>::toSecond(value); }

	ObjGrad getObjective() {return objective;}

	std::shared_ptr<std::vector<std::vector<float>>> getObjectiveTracker() {return trackedGrads;}
	std::shared_ptr<std::vector<std::vector<float>>> getStateTracker() {return trackedStates;}

	void start();

	void step();

	void finalize();

	static void enableStats() {stats = true;}

	virtual std::unique_ptr<Statistic> getStats() const = 0;


protected:

	virtual void run(push_type& ca) = 0;

	virtual Statistic* internal_getStats() = 0;

	std::shared_ptr<std::vector<std::vector<float>>> trackedStates;
	std::shared_ptr<std::vector<std::vector<float>>> trackedGrads;
	Vector state; // dof
	ObjGrad objective; // energy

	pull_type* coro;

	static bool stats;

};

}

#endif /* SOLVERBASE_H_ */
