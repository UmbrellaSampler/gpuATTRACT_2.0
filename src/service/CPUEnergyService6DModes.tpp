/*
 * CPUEnergyService6DModes.tpp
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_CPUENERGYSERVICE6DMODES_TPP_
#define SRC_CPUENERGYSERVICE6DMODES_TPP_

#include <cassert>
#include "WorkerBuffer.h"
#include "DataManager.h"
#include "DataItem.h"
#include "WorkItem.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "SimParam.h"


#include "transform_modes.h"
#include "interpolation.h"
#include "neighborlist.h"
#include "reduction_modes.h"
#include "matrixFunctions.h"
#include "RotMat.h"

// ToDo: remove
#include <iostream>

#include "CPUEnergyService6DModes.h"

namespace as {

template<typename REAL>
CPUEnergyService6DModes<REAL>::CPUEnergyService6DModes(std::shared_ptr<DataManager> dataMng) :
	CPUEnergyService<Types_6D_Modes<REAL>>(dataMng)
{}

template<typename REAL>
class CPUEnergyService6DModes<REAL>::Buffer {
public:

	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBufferRec(size_t size) {
		h_trafoLig = std::move(WorkerBuffer<REAL>(3,size));
		h_potLig = std::move(WorkerBuffer<REAL>(4,size));
	}

	void allocateBufferLig(size_t size) {
		h_trafoRec = std::move(WorkerBuffer<REAL>(3,size));
		h_defoRec = std::move(WorkerBuffer<REAL>(3,size));
		h_potRec = std::move(WorkerBuffer<REAL>(4,size));
	}

	size_t bufferSizeRec() {
		return h_trafoRec.bufferSize();
	}

	size_t bufferSizeLig() {
		return h_trafoLig.bufferSize();
	}

	WorkerBuffer<REAL> h_trafoRec;
	WorkerBuffer<REAL> h_defoRec;
	WorkerBuffer<REAL> h_trafoLig;

	WorkerBuffer<REAL> h_potRec;
	WorkerBuffer<REAL> h_potLig;
};


template<typename REAL>
auto CPUEnergyService6DModes<REAL>::createItemProcessor() -> itemProcessor_t {

	std::shared_ptr<Buffer> buffers = std::make_shared<Buffer>();

	itemProcessor_t fncObj = [this, buffers] (workItem_t* item) -> bool {
		assert(item->size() > 0);
		const auto itemSize = item->size();

		/* item pointers */
		const auto dofs = item->inputBuffer();
		const auto common = item->common();
		auto results = item->resultBuffer();

		/* get DataItem pointers */
		const auto gridRec = std::dynamic_pointer_cast<GridUnion<REAL>>(this->_dataMng->get(common->gridIdRec)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(gridRec != nullptr);

		const auto gridLig = std::dynamic_pointer_cast<GridUnion<REAL>>(this->_dataMng->get(common->gridIdLig)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(gridLig != nullptr);

		const auto lig = std::dynamic_pointer_cast<Protein<REAL>>(this->_dataMng->get(common->ligId)).get();
		assert(lig != nullptr);

		const auto rec = std::dynamic_pointer_cast<Protein<REAL>>(this->_dataMng->get(common->recId)).get();
		assert(rec != nullptr);

		const auto table = std::dynamic_pointer_cast<ParamTable<REAL>>(this->_dataMng->get(common->tableId)).get();
		assert(table != nullptr);

		const auto simParams = std::dynamic_pointer_cast<SimParam<REAL>>(this->_dataMng->get(common->paramsId)).get();
		assert(simParams != nullptr);

		if (rec->numAtoms() > buffers->bufferSizeRec()) {
			buffers->allocateBufferRec(rec->numAtoms());
		}

		if (lig->numAtoms() > buffers->bufferSizeLig()) {
			buffers->allocateBufferLig(lig->numAtoms());
		}


//		lig->print(lig->numAtoms());
//		exit(1);

		for (unsigned i = 0; i < itemSize; ++i) {
			const auto& dof = dofs[i];
			auto& enGrad = results[i];

			//invert the receptor DOF such that it points to the receptor in the ligand system
			DOF_6D_Modes<REAL> invertedRecDOF=invertDOF(dof);

			//translate the coordinates of the receptor
			rotate_translate_deform(
				rec->xPos(),
				rec->yPos(),
				rec->zPos(),
				invertedRecDOF._6D.pos,
				invertedRecDOF._6D.ang,
				rec->numAtoms(),
				rec->numModes(),
				invertedRecDOF.modesRec,
				rec->xModes(),
				rec->yModes(),
				rec->zModes(),
				buffers->h_defoRec.getX(),
				buffers->h_defoRec.getY(),
				buffers->h_defoRec.getZ(),
				buffers->h_trafoRec.getX(),//output
				buffers->h_trafoRec.getY(),
				buffers->h_trafoRec.getZ()
			); // OK

			//translate the coordinates of the Ligand
			rotate_translate(
				lig->xPos(),
				lig->yPos(),
				lig->zPos(),
				dof._6D.pos,
				dof._6D.ang,
				lig->numAtoms(),
				lig->numModes(),
				dof.modesLig,
				lig->xModes(),
				lig->yModes(),
				lig->zModes(),
				buffers->h_trafoLig.getX(),//output
				buffers->h_trafoLig.getY(),
				buffers->h_trafoLig.getZ()
			); // OK

			// Debug
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				std::cout << buffers->h_trafoLig.getX()[i] << " " << buffers->h_trafoLig.getY()[i] << " " << buffers->h_trafoLig.getZ()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);

			// calculate the forces acting on the receptor via the ligand grid in the ligand system
			potForce(
				gridLig->inner.get(),
				gridLig->outer.get(),
				rec,
				buffers->h_trafoRec.getX(),
				buffers->h_trafoRec.getY(),
				buffers->h_trafoRec.getZ(),
				buffers->h_potRec.getX(), // output
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ(),
				buffers->h_potRec.getW()
			);

			//rotate forces back into the receptor frame
			rotate_forces(invertedRecDOF._6D.ang.inv(),
				rec-> numAtoms(),
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ()
			);

			// calculate the forces acting on the ligand via the receptor grid in the receptor/global system
			potForce(
				gridRec->inner.get(),
				gridRec->outer.get(),
				lig,
				buffers->h_trafoLig.getX(),
				buffers->h_trafoLig.getY(),
				buffers->h_trafoLig.getZ(),
				buffers->h_potLig.getX(), // output
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				buffers->h_potLig.getW()
			); // OK

//			exit(EXIT_SUCCESS);
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
//			for(size_t i = 0; i < 20; ++i) {
				//std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);

			// calculate the forces acting on the receptor and the ligand in the receptor system via the neighborlist
			NLPotForce(
				gridRec->NL.get(),
				rec,
				lig,
				simParams,
				table,
				buffers->h_defoRec.getX(),
				buffers->h_defoRec.getY(),
				buffers->h_defoRec.getZ(),
				buffers->h_trafoLig.getX(),
				buffers->h_trafoLig.getY(),
				buffers->h_trafoLig.getZ(),
				buffers->h_potLig.getX(), // output
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				buffers->h_potLig.getW(),
				buffers->h_potRec.getX(), // output
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ()
			); // OK

////			// Debug
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
//			}
////			exit(EXIT_SUCCESS);




			////Reduce forces on Ligand
			PotForce_Modes<REAL> redPotForce = reducePotForce<REAL,PotForce_Modes<REAL>>(
				buffers->h_potLig.getX(),
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				buffers->h_potLig.getW(),
				lig->numAtoms()
			); // OK

//			// Debug
//			REAL x = redPotForce.pos.x;
//			REAL y = redPotForce.pos.y;
//			REAL z = redPotForce.pos.z;
//			REAL E = redPotForce.E;
//			std::cout << x << " " << y << " " << z << " " << E << std::endl;


			reduceModeForce(
				dof._6D.ang,
				buffers->h_potLig.getX(),
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				lig->xModes(),
				lig->yModes(),
				lig->zModes(),
				lig->numAtoms(),
				lig->numModes(),
				redPotForce.modesLig
				);

			correctModeForce(
				lig-> modeForce(),
				lig-> numModes(),
				redPotForce.modesLig
				);

			////Reduce forces on receptor

			reduceModeForce(
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ(),
				rec->xModes(),
				rec->yModes(),
				rec->zModes(),
				rec->numAtoms(),
				rec->numModes(),
				redPotForce.modesRec
				);

			correctModeForce(
				rec->modeForce(),
				rec-> numModes(),
				redPotForce.modesRec
				);


			//copy reduced forces
			for( int mode = 0; mode < lig->numModes(); mode++) {
				enGrad.modesLig[mode]=redPotForce.modesLig[mode];
			}
			for( int mode = 0; mode < rec->numModes(); mode++) {
				enGrad.modesRec[mode]=redPotForce.modesRec[mode];
			}
			enGrad._6D.E = redPotForce.E;
			enGrad._6D.pos = redPotForce.pos;

			enGrad._6D.ang = reduceTorque(
					lig->xPos(),
					lig->yPos(),
					lig->zPos(),
					buffers->h_potLig.getX(),
					buffers->h_potLig.getY(),
					buffers->h_potLig.getZ(),
					lig->numAtoms(),
					dof._6D.ang
			); // OK

		}

		item->setProcessed();

		return false;

	};

	return fncObj;
}

}  // namespace as




#endif /* CPUENERGYSERVICE6DMODES_TPP_ */
