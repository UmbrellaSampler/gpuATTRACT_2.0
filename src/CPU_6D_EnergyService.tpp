/*
 * CPUEnergyService.tpp
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_CPU6DENERGYSERVICE_TPP_
#define SRC_CPU6DENERGYSERVICE_TPP_

#include "CPU_6D_EnergyService.h"

#include <cassert>
#include "WorkerBuffer.h"
#include "DataManager.h"
#include "DataItem.h"
#include "WorkItem.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"
#include "SimParam.h"


#include "transform.h"
#include "interpolation.h"
#include "neighborlist.h"
#include "reduction.h"
#include "MatrixFunctions.h"
#include "RotMat.h"

// ToDo: remove
#include <iostream>

namespace as {

template<typename REAL>
class CPU_6D_EnergyService<REAL>::Buffer {
public:

	/**
	 * Allocate new Buffers with size. Old buffers are automatically deallocated;
	 */
	void allocateBuffer(size_t size) {
		h_trafoLig = std::move(WorkerBuffer<REAL>(3,size));
		h_potLig = std::move(WorkerBuffer<REAL>(4,size));
	}

	size_t bufferSize() {
		return h_trafoLig.bufferSize();
	}

	WorkerBuffer<REAL> h_trafoLig;
	WorkerBuffer<REAL> h_potLig;
};


template<typename REAL>


auto CPU_6D_EnergyService<REAL>::createItemProcessor() -> itemProcessor_t {

	std::shared_ptr<Buffer> buffers = std::make_shared<Buffer>();

	itemProcessor_t fncObj = [this,buffers] (workItem_t* item) -> bool {
		assert(item->size() > 0);
		const auto itemSize = item->size();

		/* item pointers */
		const auto dofs = item->inputBuffer();
		const auto common = item->common();
		auto results = item->resultBuffer();

		/* get DataItem pointers */
		const auto grid = std::dynamic_pointer_cast<GridUnion<REAL>>(_dataMng->get(common->gridId)).get(); // _dataMng->get() returns shared_ptr<DataItem>
		assert(grid != nullptr);

		const auto lig = std::dynamic_pointer_cast<Protein<REAL>>(_dataMng->get(common->ligId)).get();
		assert(lig != nullptr);

		const auto rec = std::dynamic_pointer_cast<Protein<REAL>>(_dataMng->get(common->recId)).get();
		assert(rec != nullptr);

		const auto table = std::dynamic_pointer_cast<ParamTable<REAL>>(_dataMng->get(common->tableId)).get();
		assert(table != nullptr);

		const auto simParams = std::dynamic_pointer_cast<SimParam<REAL>>(_dataMng->get(common->paramsId)).get();
		assert(simParams != nullptr);


		if (lig->numAtoms() > buffers->bufferSize()) {
			buffers->allocateBuffer(lig->numAtoms());
		}

//		lig->print(lig->numAtoms());
//		exit(1);

		for (unsigned i = 0; i < itemSize; ++i) {
			const auto& dof = dofs[i];
			auto& enGrad = results[i];

			rotate_translate(
					lig->xPos(),
					lig->yPos(),
					lig->zPos(),
					dof.pos,
					dof.ang,
					lig->numAtoms(),
					buffers->h_trafoLig.getX(), // output
					buffers->h_trafoLig.getY(),
					buffers->h_trafoLig.getZ()
			); // OK

			// Debug
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				std::cout << buffers->h_trafoLig.getX()[i] << " " << buffers->h_trafoLig.getY()[i] << " " << buffers->h_trafoLig.getZ()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);

			potForce(
					grid->inner.get(),
					grid->outer.get(),
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
//				std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
//			}
//			exit(EXIT_SUCCESS);

			NLPotForce(
					grid->NL.get(),
					rec,
					lig,
					simParams,
					table,
					buffers->h_trafoLig.getX(),
					buffers->h_trafoLig.getY(),
					buffers->h_trafoLig.getZ(),
					buffers->h_potLig.getX(), // output
					buffers->h_potLig.getY(),
					buffers->h_potLig.getZ(),
					buffers->h_potLig.getW()
			); // OK

////			// Debug
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
////			for(size_t i = 0; i < 20; ++i) {
//				std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
//			}
////			exit(EXIT_SUCCESS);

			PotForce<REAL> redPotForce = reducePotForce(
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

			enGrad.E = redPotForce.E;
			enGrad.pos = redPotForce.pos;

			enGrad.ang = reduceTorque(
					lig->xPos(),
					lig->yPos(),
					lig->zPos(),
					buffers->h_potLig.getX(),
					buffers->h_potLig.getY(),
					buffers->h_potLig.getZ(),
					lig->numAtoms(),
					dof.ang
			); // OK

		}

		item->setProcessed();

		return false;

	};

	return fncObj;
}

}  // namespace as


#endif /* SRC_CPU6DENERGYSERVICE_TPP_ */

// for debugging
//
//			size_t numDispl = lig->numAtoms();
//			for (size_t i = 0; i < numDispl; ++i) {
//				std::cout <<  i << " " << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
//			}
//			exit(1);
//
//			std::cout << enGrad << std::endl;
//			exit(1);
