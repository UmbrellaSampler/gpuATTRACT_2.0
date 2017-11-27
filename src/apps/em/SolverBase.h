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
#include <Eigen/Core>
#include <cassert>
#include "TypeConfiguration.h"
#include <meta.h>

namespace as {

using coro_t = boost::coroutines::coroutine<void(void)>;

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
	SolverBase() : coro(nullptr){}
	virtual ~SolverBase() { delete coro;}

	/* make object not copyable, but movealble only */
	SolverBase(const SolverBase& ) = delete;
	SolverBase& operator= (const SolverBase& ) = delete;

	SolverBase(SolverBase && rhs) {
		state = std::move(rhs.state);
		objective = std::move(rhs.objective);
		coro = std::move(rhs.coro);
		rhs.coro = nullptr;
	}

	SolverBase& operator= (SolverBase&& rhs) {
		state = std::move(rhs.state);
		objective = std::move(rhs.objective);
		coro = std::move(rhs.coro);
		rhs.coro = nullptr;
		return *this;
	}

	bool converged() {return !*coro;}
	void setState(const Vector& value) { state = value;}

//	template <typename DOFType>
//	void setState(const DOFType &value) { state = TypesConverter<DOFType, Vector>::toSecond(value);}


	template <typename DOFType>
	void setState( DOFType value) { typeConfigInput = value.getConfiguration();
	state=value.getVectorfromDOF();}


	Vector getState() {return state;}


	void setObjective(const ObjGrad& value) { objective = value; }

//	template <typename ResultType>
//	void setObjective(const ResultType& value) {	objective = TypesConverter<ResultType, ObjGrad>::toSecond(value);}

	template <typename ResultType>
	void setObjective( ResultType value) {	typeConfigResult= value.getConfiguration();
	objective=value.getVectorfromDOF();}

	ObjGrad getObjective() {return objective;}


	TypeConfiguration getInputDOFConfig(){return typeConfigInput;	}

	TypeConfiguration getResultDOFConfig(){	return typeConfigResult;}
	void start();

	void step();

	void finalize();

	static void enableStats() {stats = true;}

	virtual std::unique_ptr<Statistic> getStats() const = 0;


protected:

	virtual void run(coro_t::caller_type& ca) = 0;

	virtual Statistic* internal_getStats() = 0;

	TypeConfiguration typeConfigInput;
	TypeConfiguration typeConfigResult;

	Vector state; // dof
	ObjGrad objective; // energy

	coro_t* coro;

	static bool stats;

};
}

#endif /* SOLVERBASE_H_ */
