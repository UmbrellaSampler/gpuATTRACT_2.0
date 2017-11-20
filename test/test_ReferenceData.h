/*
 * test_ReferenceData.h
 *
 *  Created on: Nov 11, 2017
 *      Author: glenn
 */

#ifndef TEST_REFERENCEDATA_H_
#define TEST_REFERENCEDATA_H_
#include <string>
#include <iostream>

using namespace std;
//this is the data that is currently used as input for the old and new attract program . For example the DOFs, reduced pdbs for a certrain protein
const string INPUT_DATA_PATH="/home/glenn/Documents/Masterthesis/testfolder/1AVX/input/";
//In this folder the evaluation of the new units is saved. For example the evalution of error of coordinate-transformation is saved in this folder
const string EVALUATION_DATA_PATH="/home/glenn/Documents/Masterthesis/testfolder/1AVX/evaluation/";
//This folder contains data which was extracted fromt the original attract and can be used to test new units
const string TEST_DATA_PATH="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/dataset_02/";


// the following data was extracted during one cycle of transformation and reduction

//the configuration file contain the number of atoms and pivots for receptor and ligand
const string TEST_CONFIG_FILE_NAME=TEST_DATA_PATH+"config.dat";

//DOFs of ligand and receptor
const string TEST_LIG_DOF_FILE_NAME=TEST_DATA_PATH+"dofLigand.dat";
const string TEST_REC_DOF_FILE_NAME=TEST_DATA_PATH+"dofReceptor.dat";

//position of atoms before and after transformation and deformation
const string TEST_LIG_POSITION_FILE_NAME=TEST_DATA_PATH+"deformBeforeLigand.dat";
const string TEST_LIG_TRANSFORMED_POSITION_FILE_NAME=TEST_DATA_PATH+"deformAfterLigand.dat";

//modes before and after the have been read and converted by attract
const string TEST_LIG_MODE_FILE_NAME=INPUT_DATA_PATH+"modesLigand.dat";
const string TEST_LIG_CONVERTEDMODE_FILE_NAME=TEST_DATA_PATH+"writtenLigandModes.dat";

//forces acting on the ligand
const string TEST_LIG_FORCE_FILE_NAME=TEST_DATA_PATH+"forcesLigand.dat";

#endif /* TEST_REFERENCEDATA_H_ */
