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
	void allocateLigandBuffer(size_t size) {

		h_trafoLig = std::move(WorkerBuffer<REAL>(3,size));
		h_potLig = std::move(WorkerBuffer<REAL>(4,size));

	}

	void allocateReceptorBuffer(size_t size) {
		h_trafoRec = std::move(WorkerBuffer<REAL>(3,size));
		h_defoRec = std::move(WorkerBuffer<REAL>(3,size));
		h_potRec = std::move(WorkerBuffer<REAL>(4,size));
	}

	size_t bufferSizeLig() {
		return h_trafoLig.bufferSize();
	}

	size_t bufferSizeRec() {
		return h_trafoRec.bufferSize();
		}

	WorkerBuffer<REAL> h_trafoRec;
	WorkerBuffer<REAL> h_defoRec;
	WorkerBuffer<REAL> h_potRec;

	WorkerBuffer<REAL> h_trafoLig;
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


		const auto gridRec= std::dynamic_pointer_cast<GridUnion<REAL>>(this->_dataMng->get(common->gridRecId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(gridRec != nullptr);

		const auto gridLig= std::dynamic_pointer_cast<GridUnion<REAL>>(this->_dataMng->get(common->gridLigId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(gridLig != nullptr);

		const auto	lig= std::dynamic_pointer_cast<Protein<REAL>>(this->_dataMng->get(common->ligId)).get();
		assert(lig != nullptr);

		const auto rec = std::dynamic_pointer_cast<Protein<REAL>>(this->_dataMng->get(common->recId)).get();
		assert(rec != nullptr);

		const auto table = std::dynamic_pointer_cast<ParamTable<REAL>>(this->_dataMng->get(common->tableId)).get();
		assert(table != nullptr);

		const auto simParams = std::dynamic_pointer_cast<SimParam<REAL>>(this->_dataMng->get(common->paramsId)).get();
		assert(simParams != nullptr);


			if (lig->numAtoms() > buffers->bufferSizeLig()) {
				buffers->allocateLigandBuffer(lig->numAtoms());	}
			if (rec->numAtoms() > buffers->bufferSizeRec()) {
				buffers->allocateReceptorBuffer(rec->numAtoms());	}

//		lig->print(lig->numAtoms());
//		exit(1);

		for (unsigned i = 0; i < itemSize; ++i) {
			const auto& dof = dofs[i];
			auto& enGrad = results[i];
			int recIdx=0;
			int ligIdx=1;



			//invert the receptor DOF such that it points to the receptor in the ligand system
			DOF_6D_Modes<REAL> invertedRecDOF=invertDOF(dof.dof[recIdx],dof.dof[ligIdx]);

			//translate the coordinates of the receptor
			rotate_translate_deform(
				rec->xPos(),
				rec->yPos(),
				rec->zPos(),
				invertedRecDOF.pos,
				invertedRecDOF.ang,
				rec->numAtoms(),
				rec->numModes(),
				invertedRecDOF.modes,
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
			rotate_translate_deform(
				lig->xPos(),
				lig->yPos(),
				lig->zPos(),
				dof.dof[ligIdx].pos,
				dof.dof[ligIdx].ang,
				lig->numAtoms(),
				lig->numModes(),
				dof.dof[ligIdx].modes,
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
			rotate_forces(invertedRecDOF.ang.inv(),
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
			//might have to rotate some forces (cant think straight anymore, but think not because all in receptor frame which is global)
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
					buffers->h_potRec.getZ(),
					buffers->h_potRec.getW()
			); // OK

////			// Debug
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
//			}
////			exit(EXIT_SUCCESS);



			////Reduce forces on receptor

			Result_6D_Modes<REAL> redPotForceRec = reducePotForce<REAL,Result_6D_Modes<REAL>>(
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ(),
				buffers->h_potRec.getW(),
				rec->numAtoms()
			);

			reduceModeForce(
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ(),
				rec->xModes(),
				rec->yModes(),
				rec->zModes(),
				rec->numAtoms(),
				rec->numModes(),
				redPotForceRec.modes
				);

			correctModeForce(
				rec->modeForce(),
				rec-> numModes(),
				redPotForceRec.modes
				);


			redPotForceRec.ang = reduceTorque(
				buffers->h_defoRec.getX(),
				buffers->h_defoRec.getY(),
				buffers->h_defoRec.getZ(),
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ(),
				rec->numAtoms(),
				dof.dof[recIdx].ang
			); // OK

			enGrad.dof.push_back(redPotForceRec);

			////Reduce forces on ligand

			Result_6D_Modes<REAL> redPotForceLig = reducePotForce<REAL,Result_6D_Modes<REAL>>(
				buffers->h_potLig.getX(),
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				buffers->h_potLig.getW(),
				lig->numAtoms()
			); // OK

			reduceModeForce(
				dof.dof[ligIdx].ang,
				buffers->h_potLig.getX(),
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				lig->xModes(),
				lig->yModes(),
				lig->zModes(),
				lig->numAtoms(),
				lig->numModes(),
				redPotForceLig.modes
				);

			correctModeForce(
				lig-> modeForce(),
				lig-> numModes(),
				redPotForceLig.modes
				);


			redPotForceLig.ang = reduceTorque(
				lig->xPos(),
				lig->yPos(),
				lig->zPos(),
				buffers->h_potLig.getX(),
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				lig->numAtoms(),
				dof.dof[ligIdx].ang
			); // OK

			// doing it this way it is important which DOF is pushed back first!!!
			enGrad.dof.push_back(redPotForceLig);

		}

		item->setProcessed();

		return false;

	};

	return fncObj;
}

}  // namespace as




#endif /* CPUENERGYSERVICE6DMODES_TPP_ */
