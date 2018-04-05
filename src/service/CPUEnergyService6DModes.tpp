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
#include "neighborlist_modes.h"
#include "reduction_modes.h"
#include "matrixFunctions.h"
#include "RotMat.h"

// ToDo: remove
#include <iostream>
#include "readFile.h"
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
	void allocateBufferLig(size_t size) {
		h_defoLig = std::move(WorkerBuffer<REAL>(3,size));
		h_trafoLig = std::move(WorkerBuffer<REAL>(3,size));
		h_potLig = std::move(WorkerBuffer<REAL>(4,size));
	}

	void allocateBufferRec(size_t size) {
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
	WorkerBuffer<REAL> h_defoLig;

	WorkerBuffer<REAL> h_potRec;
	WorkerBuffer<REAL> h_potLig;
};


template <typename REAL>
Vec3<REAL> getDifferenceOfVector(Vec3<REAL> ligPos, Vec3<REAL> recPos){
	return ligPos - recPos;
}

template <typename REAL>
Vec3<REAL> normedRatioOfVector(Vec3<REAL> vec_denominator, Vec3<REAL> vec_enumerator){
	Vec3<REAL> ratio
	{
	abs(abs( vec_enumerator.x / vec_denominator.x ) - 1.0 ),
	abs(abs( vec_enumerator.y / vec_denominator.y ) - 1.0 ),
	abs(abs( vec_enumerator.z / vec_denominator.z ) - 1.0 )
	};
	return ratio;
}

template <typename REAL>
bool ifREALsmallerVector(Vec3<REAL> vec, REAL epsilon){
	if( abs(vec.x ) > epsilon ||
		abs(vec.y ) > epsilon ||
		abs(vec.z ) > epsilon )
	{
		return true;
	}
	else
	{
		return false;
	}
}
template < typename REAL>
void compareTwoVec3(Vec3<REAL> vec1, Vec3<REAL> vec2, REAL epsilon)
{
	Vec3<REAL> normedRatio = normedRatioOfVector(vec1, vec2);
	if( ifREALsmallerVector(normedRatio, epsilon) )
	{
		std::cout << getDifferenceOfVector(vec1, vec2) << std::endl;
		}
}

template <typename T1, typename T2, typename REAL>
void compareTwo3colArrays(T1 const& array1, T2 const& array2, unsigned size_row, REAL epsilon){
	for( unsigned row; row < size_row; row++)
	{
		Vec3<REAL> vec1{array1[row][0],
						array1[row][1],
			            array1[row][2]};
		Vec3<REAL> vec2{array2[row][0],
						array2[row][1],
						array2[row][2]};
		compareTwoVec3(vec1, vec2, epsilon);
	}
}


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

			Vec3<REAL> pr(32.823001861572266, -6.3949999809265137, 23.483999252319336);
			Vec3<REAL> pl(51.143001556396484, 7.6799998283386230, 38.110000610351562);
			DOF_6D_Modes<REAL> test = dof;
			test._6D.pos = test._6D.pos + pr - pl;


			//std::cout << test << std::endl;

			//invert the receptor DOF such that it points to the receptor in the ligand system
			DOF_6D_Modes<REAL> invertedRecDOF=invertDOF(dof);

//			Vec3<REAL> rec_ang(0.0);
//			Vec3<REAL> rec_pos(0.0);
//	 		rotMatr = euler2rotmat(rec_ang.x, rec_ang.y, rec_ang.z);
//	 		rotMatl = euler2rotmat(ang.x, ang.y, ang.z);
//	 		rotMatrInv  = rotMatr.getInv();
//	 		rotMatlInv  = rotMatl.getInv();
//	 		rotMatd 	= rotMatl * rotMatrInv;
//	 		Vec3<REAL> trans_d = dof._6D.pos - rec_pos;
//	 		trans_d = rotMatd * trans_d;
//	 		Vec3<REAL> pivot_d = pl - pr;
//	 		pivot_d = rotMatd * pivot_d;
//	 		pivot_d = pivot_d + pr;
//	 		Vec3<REAL> trans = pivot_d + trans_d;


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
				buffers->h_trafoRec.getZ(),
				0
			); // OK


			rotate_translate_deform(
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
				buffers->h_defoLig.getX(),
				buffers->h_defoLig.getY(),
				buffers->h_defoLig.getZ(),
				buffers->h_trafoLig.getX(),//output
				buffers->h_trafoLig.getY(),
				buffers->h_trafoLig.getZ(),
				1
			); // OK



//			std::vector<std::vector<REAL>> ligandPivoOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/deformBeforeLigandPivotized_0000.dat" );
//			std::vector<std::vector<REAL>> receptorPivoOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/deformBeforeReceptorPivotized_0000.dat" );
//
//			std::vector<std::vector<REAL>> receptorModesOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/writtenReceptorModes.dat" );
//			std::vector<std::vector<REAL>> ligandModesOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/writtenLigandModes.dat" );
//
//			std::vector<std::vector<REAL>> ligandTransOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/deformAfterLigand_0000.dat" );
//			std::vector<std::vector<REAL>> receptorTransOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/deformAfterReceptor_0000.dat" );
//
//			std::vector<std::vector<REAL>> ligandDefoOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/deformBeforeLigand_0000.dat" );
//			std::vector<std::vector<REAL>> receptorDefoOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/deformBeforeReceptor_0000.dat" );
//
//			std::vector<std::vector<REAL>> ligandDefoAmplitudeOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/lig_deform_amplitude.dat" );
//			std::vector<std::vector<REAL>> receptorDefoAmplitudeOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/rec_deform_amplitude.dat" );
//
//			Vec3<REAL> pr(32.823001861572266, -6.3949999809265137, 23.483999252319336);
//			Vec3<REAL> pl(51.143001556396484, 7.6799998283386230, 38.110000610351562);
//
//			REAL epsmod = 0.001;
//			unsigned modeIdx = 0;
//			unsigned realIdx = 0;
//			std::cout << "ligand sizes. size orig" << ligandTransOrig.size () << "neew size" << lig->numAtoms() << std::endl;
//			std::cout << "receptor sizes. orig" << receptorTransOrig.size () << "neew size" << rec->numAtoms() << std::endl;
//
//
//			std::cout << "\n Compare ligand Modes" << std::endl;
//			for(size_t mode = 0; mode < 1; ++mode)
//			{
//				for(size_t i = 0; i < lig->numAtoms(); ++i)
//				{
//					realIdx = (lig->numAtoms() + 1) * mode + i + 1;
//					Vec3<REAL> origModeLig(ligandModesOrig[realIdx][0], ligandModesOrig[realIdx][1], ligandModesOrig[realIdx][2] );
//					Vec3<REAL> newModeLig(lig->xModes()[i * 5 + mode],lig->yModes()[i * 5 + mode], lig->zModes()[i * 5 + mode] );
//					compareTwoVec3(origModeLig, newModeLig, epsmod);
//				}
//			}
//
//
//			std::cout << "\n Compare receptor modes" << std::endl;
//			for(size_t mode = 0; mode < 5; ++mode)
//			{
//				for(size_t i = 0; i < rec->numAtoms(); ++i)
//				{
//					realIdx = (rec->numAtoms() + 1) * mode + i + 1;
//				//	Vec3<REAL> modeRec(receptorModesOrig[realIdx][0], receptorModesOrig[realIdx][1], receptorModesOrig[realIdx][2] );
//					Vec3<REAL> origModeRec(receptorModesOrig[realIdx][0],receptorModesOrig[realIdx][1], receptorModesOrig[realIdx][2] );
//					Vec3<REAL> newModeRec(rec->xModes()[i * 5 + mode],rec->yModes()[i * 5 + mode], rec->zModes()[i * 5 + mode] );
//					compareTwoVec3(origModeRec, newModeRec, epsmod);
//				}
//			}
//
//			REAL epsloc = 0.005;
//
//			std::cout << "\n Compare receptor Location" << std::endl;
//			for(size_t i = 0; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> diffNewRec( buffers->h_trafoRec.getX()[i]- buffers->h_defoLig.getX()[0],
//									   buffers->h_trafoRec.getY()[i]- buffers->h_defoLig.getY()[0],
//									   buffers->h_trafoRec.getZ()[i]- buffers->h_defoLig.getZ()[0] );
//				Vec3<REAL> diffOrigRec(receptorTransOrig[i][0] - ligandDefoOrig[0][0],
//									   receptorTransOrig[i][1] - ligandDefoOrig[0][1],
//									   receptorTransOrig[i][2] - ligandDefoOrig[0][2]);
//				compareTwoVec3(diffOrigRec, diffNewRec, epsloc);
//			}
//
//			std::cout << "\n Compare ligand Location" << std::endl;
//			for(size_t i = 0; i < lig->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> diffNewLig ( buffers->h_trafoLig.getX()[i] - buffers->h_defoRec.getX()[10] + pr.x,
//										buffers->h_trafoLig.getY()[i] - buffers->h_defoRec.getY()[10] + pr.y,
//										buffers->h_trafoLig.getZ()[i] - buffers->h_defoRec.getZ()[10] + pr.z);
//				//diffNewLig = diffNewLig - pr;
//				Vec3<REAL> diffOrigLig( ligandTransOrig[i][0] - receptorPivoOrig[10][0],
//										ligandTransOrig[i][1] - receptorPivoOrig[10][1],
//										ligandTransOrig[i][2] - receptorPivoOrig[10][2]);
//
//				compareTwoVec3(diffOrigLig, diffNewLig, epsloc);
//			}
//
//
//			std::cout << "\n Compare receptor internal location Difference deformation" << std::endl;
//			for(size_t i = 1; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> diffNewRec( buffers->h_defoRec.getX()[i] - buffers->h_defoRec.getX()[0],
//									   buffers->h_defoRec.getY()[i] - buffers->h_defoRec.getY()[0],
//									   buffers->h_defoRec.getZ()[i] - buffers->h_defoRec.getZ()[0] );
//				Vec3<REAL> diffOrigRec(receptorDefoOrig[i][0] - receptorDefoOrig[0][0],
//									   receptorDefoOrig[i][1] - receptorDefoOrig[0][1],
//									   receptorDefoOrig[i][2] - receptorDefoOrig[0][2]);
//				//diffNewRec = diffNewRec - pl;
//				if( abs(abs( diffNewRec.x / diffOrigRec.x ) - 1.0 ) > epsloc ||
//					abs(abs( diffNewRec.y / diffOrigRec.y ) - 1.0 ) > epsloc  ||
//					abs(abs( diffNewRec.z / diffOrigRec.z ) - 1.0 ) > epsloc )
//					compareTwoVec3(diffOrigRec, diffNewRec, epsloc);
//			}
//
//			std::cout << "\n Compare receptor internal location Difference transformation" << std::endl;
//			for(size_t i = 1; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> diffNewRec( buffers->h_trafoRec.getX()[i] - buffers->h_trafoRec.getX()[0],
//									   buffers->h_trafoRec.getY()[i] - buffers->h_trafoRec.getY()[0],
//									   buffers->h_trafoRec.getZ()[i] - buffers->h_trafoRec.getZ()[0] );
//				Vec3<REAL> diffOrigRec(receptorTransOrig[i][0] - receptorTransOrig[0][0],
//									   receptorTransOrig[i][1] - receptorTransOrig[0][1],
//									   receptorTransOrig[i][2] - receptorTransOrig[0][2]);
//				compareTwoVec3(diffOrigRec, diffNewRec, epsloc);
//			}
//
//
//
//
//			std::cout << "\n Compare ligand internal location Difference" << std::endl;
//			for(size_t i = 1; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> diffNewLig( buffers->h_trafoLig.getX()[i] - buffers->h_trafoLig.getX()[0],
//									   buffers->h_trafoLig.getY()[i] - buffers->h_trafoLig.getY()[0],
//									   buffers->h_trafoLig.getZ()[i] - buffers->h_trafoLig.getZ()[0] );
//				Vec3<REAL> diffOrigLig(ligandTransOrig[i][0] - ligandTransOrig[0][0],
//									   ligandTransOrig[i][1] - ligandTransOrig[0][1],
//									   ligandTransOrig[i][2] - ligandTransOrig[0][2]);
//				//diffNewRec = diffNewRec - pl;
//				compareTwoVec3(diffOrigLig, diffNewLig, epsloc);
//			}
//
//
//			REAL epslocpivo = 0.01;
//
//
//
//			std::cout << "\n Compare ligand original pivotized coordinates" << std::endl;
//			for(size_t i = 0; i < lig->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> pivoNewLig ( lig->xPos()[i],
//										lig->yPos()[i],
//										lig->zPos()[i]);
//
//				Vec3<REAL> pivoOrigLig( ligandPivoOrig[i][0],
//										ligandPivoOrig[i][1],
//										ligandPivoOrig[i][2]);
//				compareTwoVec3(pivoOrigLig, pivoNewLig, epslocpivo);
//			}
//
//			std::cout << "\n Compare receptor original pivotized coordinates" << std::endl;
//			for(size_t i = 0; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> pivoNewRec ( rec->xPos()[i],
//										rec->yPos()[i],
//										rec->zPos()[i]);
//
//				Vec3<REAL> pivoOrigRec( receptorPivoOrig[i][0],
//										receptorPivoOrig[i][1],
//										receptorPivoOrig[i][2]);
//				if( 	abs(abs( pivoNewRec.x / pivoOrigRec.x ) - 1.0 ) > epsloc ||
//						abs(abs( pivoNewRec.y / pivoOrigRec.y ) - 1.0 ) > epsloc ||
//						abs(abs( pivoNewRec.z / pivoOrigRec.z ) - 1.0 ) > epsloc )
//				{
////					std::cout << "gpuattract pivo " << i << " " << pivoNewRec << std::endl;
////					std::cout << "original   pivo " << i << " " << pivoOrigRec << std::endl;
////					std::cout <<  std::endl;
//				}
//			}
//
//
//
//			std::cout << "\n Compare ligand deformation amplitude" << std::endl;
//			for(size_t i = 0; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> defoAmplNewLig (buffers->h_defoLig.getX()[i],
//										   	buffers->h_defoLig.getY()[i],
//										   	buffers->h_defoLig.getZ()[i]);
//
//				Vec3<REAL> defoAmplOrigLig(   ligandDefoAmplitudeOrig[i][0],
//										   	   ligandDefoAmplitudeOrig[i][1],
//										   	   ligandDefoAmplitudeOrig[i][2]);
//				if( 	abs(abs( defoAmplNewLig.x / defoAmplOrigLig.x ) - 1.0 ) > 0.1 ||
//						abs(abs( defoAmplNewLig.y / defoAmplOrigLig.y ) - 1.0 ) > 0.1 ||
//						abs(abs( defoAmplNewLig.z / defoAmplOrigLig.z ) - 1.0 ) > 0.1 )
//				{
//					std::cout << "gpuattract deformation ampplitude " << i << " " << defoAmplNewLig << std::endl;
//					std::cout << "original   deformation ampplitude " << i << " " << defoAmplOrigLig << std::endl;
//					std::cout <<  std::endl;
//				}
//			}
//
//			std::cout << "\n Compare receptor deformation amplitude" << std::endl;
//			for(size_t i = 0; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> defoAmplNewRec ( 	buffers->h_defoRec.getX()[i],
//												buffers->h_defoRec.getY()[i],
//												buffers->h_defoRec.getZ()[i]);
//
//				Vec3<REAL> defoAmplOrigRec( 	receptorDefoAmplitudeOrig[i][0],
//												receptorDefoAmplitudeOrig[i][1],
//												receptorDefoAmplitudeOrig[i][2]);
//				if( 	abs(abs( defoAmplNewRec.x / defoAmplOrigRec.x ) - 1.0 ) > 0.1 ||
//						abs(abs( defoAmplNewRec.y / defoAmplOrigRec.y ) - 1.0 ) > 0.1 ||
//						abs(abs( defoAmplNewRec.z / defoAmplOrigRec.z ) - 1.0 ) > 0.1 )
//				{
////					std::cout << "gpuattract deformation amplitude " << i << " " << defoAmplNewRec << std::endl;
////					std::cout << "original   deformation amplitude " << i << " " << defoAmplOrigRec << std::endl;
////					std::cout << abs( defoAmplNewRec.x / defoAmplOrigRec.x ) << " " << abs( defoAmplNewRec.y / defoAmplOrigRec.y )<< " " <<abs( defoAmplNewRec.z / defoAmplOrigRec.z )<< std::endl;
//				}
//			}
//
//
//			std::cout << "\n Compare ligand deformed coordinates" << std::endl;
//			for(size_t i = 0; i < lig->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> defoNewLig ( buffers->h_defoLig.getX()[i],
//										buffers->h_defoLig.getY()[i],
//										buffers->h_defoLig.getZ()[i]);
//				defoNewLig = defoNewLig + pl;
//				Vec3<REAL> defoOrigLig( ligandDefoOrig[i][0],
//										ligandDefoOrig[i][1],
//										ligandDefoOrig[i][2]);
//				if( 	abs(abs( defoNewLig.x / defoOrigLig.x ) - 1.0 ) > epsloc ||
//						abs(abs( defoNewLig.y / defoOrigLig.y ) - 1.0 ) > epsloc ||
//						abs(abs( defoNewLig.z / defoOrigLig.z ) - 1.0 ) > epsloc )
//				{
////					std::cout << "gpuattract defo  " << i << " " << defoNewLig << std::endl;
////					std::cout << "original defo " << i << " " << defoOrigLig << std::endl;
////					std::cout <<  std::endl;
//				}
//			}
//
//
//			std::cout << "\n Compare receptor deformed coordinates" << std::endl;
//			for(size_t i = 0; i < rec->numAtoms(); ++i)
//			{
//
//				Vec3<REAL> defoNewRec ( buffers->h_defoRec.getX()[i],
//										buffers->h_defoRec.getY()[i],
//										buffers->h_defoRec.getZ()[i]);
//				defoNewRec = defoNewRec + pr;
//				Vec3<REAL> defoOrigRec( receptorDefoOrig[i][0],
//										receptorDefoOrig[i][1],
//										receptorDefoOrig[i][2]);
//			if( abs(abs( defoNewRec.x / defoOrigRec.x ) - 1.0 ) > epsloc ||
//				abs(abs( defoNewRec.y / defoOrigRec.y ) - 1.0 ) > epsloc ||
//				abs(abs( defoNewRec.z / defoOrigRec.z ) - 1.0 ) > epsloc )
//				{
////					std::cout << "gpuattract defo  " << i << " " << defoNewRec << std::endl;
////					std::cout << "original defo " << i << " " << defoOrigRec << std::endl;
////					std::cout <<  std::endl;
//				}
//			}





			//translate the coordinates of the Ligand
//			rotate_translate(
//				lig->xPos(),
//				lig->yPos(),
//				lig->zPos(),
//				dof._6D.pos,
//				dof._6D.ang,
//				lig->numAtoms(),
//				lig->numModes(),
//				dof.modesLig,
//				lig->xModes(),
//				lig->yModes(),
//				lig->zModes(),
//				buffers->h_trafoLig.getX(),//output
//				buffers->h_trafoLig.getY(),
//				buffers->h_trafoLig.getZ()
//			); // OK

//			// Debug
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

			rotate_forces(invertedRecDOF._6D.ang,
				rec-> numAtoms(),
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ()
			);
//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
//			for(size_t i = 0; i < 20; ++i) {
//				std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << " " << buffers->h_potLig.getW()[i] << std::endl;
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


			//			for(size_t i = 0; i < lig->numAtoms(); ++i) {
			for(size_t i = 0; i < 20; ++i) {
				//	std::cout << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << std::endl;//<< " " << buffers->h_potLig.getW()[i] << std::endl;
			}

//			NLPotForce(
//				gridRec->NL.get(),
//				rec,
//				lig,
//				simParams,
//				table,
//				buffers->h_defoRec.getX(),
//				buffers->h_defoRec.getY(),
//				buffers->h_defoRec.getZ(),
//				buffers->h_trafoLig.getX(),
//				buffers->h_trafoLig.getY(),
//				buffers->h_trafoLig.getZ(),
//				buffers->h_potLig.getX(), // output
//				buffers->h_potLig.getY(),
//				buffers->h_potLig.getZ(),
//				buffers->h_potLig.getW()
//				);
//
//			NLPotForce(
//				gridLig->NL.get(),
//				lig,
//				rec,
//				simParams,
//				table,
//				buffers->h_defoLig.getX(),
//				buffers->h_defoLig.getY(),
//				buffers->h_defoLig.getZ(),
//				buffers->h_trafoRec.getX(),
//				buffers->h_trafoRec.getY(),
//				buffers->h_trafoRec.getZ(),
//				buffers->h_potRec.getX(), // output
//				buffers->h_potRec.getY(),
//				buffers->h_potRec.getZ(),
//				buffers->h_potRec.getW()
//				);




//			NLPotForce(
//				gridLig->NL.get(),
//				lig,
//				rec,
//				simParams,
//				table,
//				buffers->h_defoLig.getX(),
//				buffers->h_defoLig.getY(),
//				buffers->h_defoLig.getZ(),
//				buffers->h_trafoRec.getX(),
//				buffers->h_trafoRec.getY(),
//				buffers->h_trafoRec.getZ(),
//				buffers->h_potRec.getX(), // output
//				buffers->h_potRec.getY(),
//				buffers->h_potRec.getZ(),
//				buffers->h_potRec.getW(),
//				buffers->h_potLig.getX(), // output
//				buffers->h_potLig.getY(),
//				buffers->h_potLig.getZ()
//			); // OK

////			// Debug

//			std::vector<std::vector<REAL>> ligandForcesOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/forcesLigand_0000.dat" );
//			std::vector<std::vector<REAL>> receptorForcesOrig = readArray<REAL>( "/home/glenn/Documents/Masterthesis/testfolder/Uwestestfile/attract_2018_02_23/orig_attract_data/forcesReceptor_0000.dat" );
//
//
//			REAL eps = 0.01;
//
//			std::cout << " ligand sizes. size orig" << ligandForcesOrig.size () << "neew size" << lig->numAtoms() << std::endl;
//			std::cout << " receptor sizes. orig" << receptorForcesOrig.size () << "neew size" << rec->numAtoms() << std::endl;
//			Vec3<REAL> diffSumRec(0.0);
//			Vec3<REAL> diffSumLig(0.0);
//
//			assert( lig->numAtoms() ==   ligandForcesOrig.size() );
//			assert( rec->numAtoms() == receptorForcesOrig.size() );
//			std::cout << "\n Compare ligand Forces" << std::endl;
//
//			for(size_t i = 0; i < lig->numAtoms(); ++i)
//			{
//				if( abs( abs( buffers->h_potLig.getX()[i] / ligandForcesOrig[i][0]) -1.0 ) > eps || abs(abs( buffers->h_potLig.getY()[i] / ligandForcesOrig[i][1]) -1.0 ) > eps || abs(abs( buffers->h_potLig.getZ()[i] / ligandForcesOrig[i][2]) - 1.0 ) > eps )
//				{
//					//std::cout << "gpuattractForces  " << i << " " << buffers->h_potLig.getX()[i] << " " << buffers->h_potLig.getY()[i] << " " << buffers->h_potLig.getZ()[i] << std::endl;
//					//std::cout << "original Forces   " << i << " " << ligandForcesOrig[i][0] << " " << ligandForcesOrig[i][1] << " " << ligandForcesOrig[i][2] << std::endl;
//
//					diffSumLig.x += buffers->h_potLig.getX()[i] - ligandForcesOrig[i][0];
//					diffSumLig.y += buffers->h_potLig.getY()[i] - ligandForcesOrig[i][1];
//					diffSumLig.z += buffers->h_potLig.getZ()[i] - ligandForcesOrig[i][2];
//					//std::cout << "difference     " << i << " " << diffSumLig << std::endl;
//					//std::cout <<  std::endl;
//
//				}
//			}
//
//
//
//			std::cout << "\n Compare receptor Forces" << std::endl;
//			for(size_t i = 0; i < rec->numAtoms(); ++i)
//			{
//				if( abs(abs( buffers->h_potRec.getX()[i] / receptorForcesOrig[i][0]) - 1.0 ) > eps ||
//					abs(abs( buffers->h_potRec.getY()[i] / receptorForcesOrig[i][1]) -1.0 ) > eps ||
//					abs(abs( buffers->h_potRec.getZ()[i] / receptorForcesOrig[i][2]) -1.0 ) > eps )
//				{
//				//	std::cout << "gpuattractForces " << i << " " << buffers->h_potRec.getX()[i] << " " << buffers->h_potRec.getY()[i] << " " << buffers->h_potRec.getZ()[i] << std::endl;
//				//	std::cout << "original Forces  " << i << " " << receptorForcesOrig[i][0] << " " << receptorForcesOrig[i][1] << " " << receptorForcesOrig[i][2] << std::endl;
//
//
//					diffSumRec.x = buffers->h_potRec.getX()[i] - receptorForcesOrig[i][0];
//					diffSumRec.y = buffers->h_potRec.getY()[i] - receptorForcesOrig[i][1];
//					diffSumRec.z = buffers->h_potRec.getZ()[i] - receptorForcesOrig[i][2];
//				//	std::cout << "difference     " << i << " " << diffSumRec << std::endl;
//				//	std::cout <<  std::endl;
//				}
//			}
//
//			std::cout << "forces difference receptor" << diffSumRec << std::endl;
//			std::cout << "forces difference ligand " << diffSumLig << std::endl;
			//exit(EXIT_SUCCESS);




			////Reduce forces on Ligand
			PotForce_Modes<REAL> redPotForce = reducePotForce<REAL,PotForce_Modes<REAL>>(
				buffers->h_potLig.getX(),
				buffers->h_potLig.getY(),
				buffers->h_potLig.getZ(),
				buffers->h_potLig.getW(),
				lig->numAtoms()
			); // OK


			PotForce_Modes<REAL> redPotForceReceptor = reducePotForce<REAL,PotForce_Modes<REAL>>(
				buffers->h_potRec.getX(),
				buffers->h_potRec.getY(),
				buffers->h_potRec.getZ(),
				buffers->h_potRec.getW(),
				rec->numAtoms()
			); // OK

//			redPotForce.pos += redPotForceReceptor.pos;
//			// Debug
//			REAL x = redPotForce.pos.x;
//			REAL y = redPotForce.pos.y;
//			REAL z = redPotForce.pos.z;
//			REAL E = redPotForce.E;
//			std::cout << x << " " << y << " " << z << " " << E << std::endl;


			reduceModeForce(
				//invertedRecDOF._6D.ang,
				Vec3<REAL> (0.0),
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
			// std::cout << " \n\n";
           // std::cout << " lig mode ";
			correctModeForce(
				lig-> modeForce(),
				lig-> numModes(),
				dof.modesLig,
				redPotForce.modesLig
				);

			////Reduce forces on receptor
			//std::cout << " rec mode ";
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
				dof.modesRec,
				redPotForce.modesRec
				);


			//copy reduced forces
			//std::cout << "print out mode force ligand"<< std::endl;
			for( int mode = 0; mode < lig->numModes(); mode++) {
				enGrad.modesLig[mode]=redPotForce.modesLig[mode];
				//enGrad.modesLig[mode]=0;
				//	std::cout << redPotForce.modesLig[mode] << " ";
			}
			//std::cout << "print out mode force receptor"<< std::endl;

			for( int mode = 0; mode < rec->numModes(); mode++) {
				enGrad.modesRec[mode]=redPotForce.modesRec[mode];
				//	std::cout << redPotForce.modesRec[mode] << " ";
				//enGrad.modesRec[mode]=0;
			}
			//std::cout <<"ligand mode energy ";
			double modeEnergyLigand = getModeEngergy(lig->modeForce(),
					lig->numModes(),
					dof.modesLig
					);
			//std::cout <<"receptor mode energy ";
			double modeEnergyReceptor = getModeEngergy(rec->modeForce(),
					rec->numModes(),
					dof.modesRec
					);

			enGrad._6D.E = redPotForce.E + modeEnergyReceptor + modeEnergyLigand;
			//std::cout << "the total energy is: "<< redPotForce.E << " the energy of the lig: "<< modeEnergyLigand << " E rec: " << modeEnergyReceptor<< std::endl;
			enGrad._6D.pos = redPotForce.pos ;//- redPotForceReceptor.pos;

			//std::cout << " print out forconstanst";
			for( int i = 0; i < 5; i++){
			//	std::cout << "i: " << i << "lig : "<< lig->modeForce()[i] << "rec: " << rec->modeForce()[i]<< std::endl;
			}
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


			for( int i = 0; i< 10; i++){
				//std::cout << buffers->h_potRec.getX()[i] << " " << buffers->h_potRec.getY()[i] << " " << buffers->h_potRec.getZ()[i]<< std::endl;

			}
			Vec3<REAL> angForce = reduceTorque(
					rec->xPos(),
					rec->yPos(),
					rec->zPos(),
					buffers->h_potRec.getX(),
					buffers->h_potRec.getY(),
					buffers->h_potRec.getZ(),
					rec->numAtoms(),
					Vec3<REAL>(0.0)
						); // OK
//			std::cout << std::endl;
//			std::cout << "Ligand Deltas  " << std::endl;
//
//			std::cout << "ROTATIONAL     " << enGrad._6D.ang << std::endl;
//			std::cout << "TRANSLATIONAL  " << redPotForce.pos << std::endl;
//			std::cout << "MODES Lig         " ;
//			for( int mode = 0; mode < lig->numModes(); mode++) {
//				std::cout << redPotForce.modesLig[mode] << " ";
//			}
//			std::cout << std::endl;
//			std::cout << std::endl;
//			std::cout << "Receptor Deltas" << std::endl;

		//	std::cout << "ROTATIONAL     " << angForce << std::endl;
		//	std::cout << "TRANSLATIONAL  " << redPotForceReceptor.pos << std::endl;
//			std::cout << "MODES Rec         " ;
//			for( int mode = 0; mode < rec->numModes(); mode++) {
//				std::cout << redPotForce.modesRec[mode] << " ";
//			}
//			std::cout << std::endl;
//			double what = 1;
		}

		item->setProcessed();

		return false;

	};

	return fncObj;
}

}  // namespace as




#endif /* CPUENERGYSERVICE6DMODES_TPP_ */
