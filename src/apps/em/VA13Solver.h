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

#ifndef VA13SOLVER_H_
#define VA13SOLVER_H_

#include "SolverBase.h"

namespace as {

struct VA13Statistic : public Statistic {

	virtual Statistic* getCopy() const override {
		return static_cast<Statistic*> (new VA13Statistic(*this));
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

class VA13Solver : public SolverBase {
public:
	VA13Solver() : SolverBase() {}
	virtual ~VA13Solver() {};

	VA13Solver(const VA13Solver& ) = delete;
	VA13Solver& operator= (const VA13Solver& ) = delete;

	VA13Solver(VA13Solver &&) = default;
	VA13Solver& operator= (VA13Solver&& ) = default;

	std::unique_ptr<Statistic> getStats() const override {
		return std::unique_ptr<Statistic> (statistic.getCopy());
	}

	struct Options {
		/* Solver Options */
		unsigned maxFunEval = 500;
		unsigned minimizeRotation = 1;  //switch on(1)/off(!1) minimization of rotational degrees of freedom
		unsigned minimizeTranslation = 1; //switch on(1)/off(!1) minimization of translational degrees of freedom
		unsigned minimizeModes = 1; //switch on(1)/off(!1) minimization of mode degrees of freedom
		bool trackStates = true; //switch on/off tracking of states. which are appended to trackedStates
		bool trackGradients = true; //switch on/off tracking of gradients. which are appended to trackedGrads
	};

	static void setOptions(Options opt) {settings = opt;}

	class FortranSmuggler {
	public:
		FortranSmuggler (push_type& _coro, Vector& _state, ObjGrad& _objective,std::shared_ptr<std::vector<std::vector<float>>> _trackedStates,
				std::shared_ptr<std::vector<std::vector<float>>> _trackedGrads,Options opt ):
			coro(_coro),
			state(_state),
			objective(_objective),
			trackedStates(_trackedStates),
			trackedGrads(_trackedGrads),
			options(opt)
		{}

		Vector& state_ref() { return state; }
		ObjGrad& objective_ref() { return objective; }
		void push_state(std::vector<float> state) { trackedStates->push_back(state);}
		void push_grad(std::vector<float> grad) { trackedGrads->push_back(grad); }
		void call_coro() {
			coro();
		};
		Options getOptions(){return options;}

	private:
		push_type& coro;
		Vector& state;
		ObjGrad& objective;
		Options options;
		std::shared_ptr<std::vector<std::vector<float>>> trackedStates;
		std::shared_ptr<std::vector<std::vector<float>>> trackedGrads;
	};


private:

	void run(push_type& ca) override;

	/* solver options */
	static Options settings;

	/* Statistics */
	VA13Statistic statistic;

	virtual Statistic* internal_getStats() override {
		return static_cast<Statistic*>(&statistic);
	}

};

}


#endif /* VA13SOLVER_H_ */
