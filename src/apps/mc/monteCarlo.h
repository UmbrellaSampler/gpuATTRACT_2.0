/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2018 Glenn Glashagen
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

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "SolverBase.h"

namespace as {

struct MCStatistic : public Statistic {

	virtual Statistic* getCopy() const override {
		return static_cast<Statistic*> (new MCStatistic(*this));
	}

	virtual std::ostream& print(std::ostream& stream) const override {
		using namespace std;
		int precisionSetting = stream.precision( );
		ios::fmtflags flagSettings = stream.flags();
		stream.setf(ios::scientific);
		stream.precision(5);

		stream << numRequests << endl;

		stream.precision(precisionSetting);
		stream.flags(flagSettings);
		return stream;
	}
};

class MCSolver : public SolverBase {
public:
	MCSolver() : SolverBase() {}
	virtual ~MCSolver() {};

	MCSolver(const MCSolver& ) = delete;
	MCSolver& operator= (const MCSolver& ) = delete;

	MCSolver(MCSolver &&) = default;
	MCSolver& operator= (MCSolver&& ) = default;

	std::unique_ptr<Statistic> getStats() const override {
		return std::unique_ptr<Statistic> (statistic.getCopy());
	}

	struct Options {
		/* Solver Options */
		unsigned maxFunEval = 500;
	};

	static void setOptions(Options opt) {settings = opt;}


	class FortranSmuggler {
	public:
		FortranSmuggler (push_type& _coro, Vector& _state, ObjGrad& _objective):
			coro(_coro),
			state(_state),
			objective(_objective)
		{}

		Vector& state_ref() { return state; }
		ObjGrad& objective_ref() { return objective; }
		void call_coro() {
			coro();
		};

	private:
		push_type& coro;
		Vector& state;
		ObjGrad& objective;
	};


private:

	void run(push_type& ca) override;

	/* solver options */
	static Options settings;

	/* Statistics */
	MCStatistic statistic;

	virtual Statistic* internal_getStats() override {
		return static_cast<Statistic*>(&statistic);
	}

};

}


#endif /* MCSOLVER_H_ */
