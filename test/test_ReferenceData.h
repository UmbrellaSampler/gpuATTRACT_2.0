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

const string EVALUATION_DATA_PATH="/home/glenn/Documents/Masterthesis/testfolder/1AVX/evaluation/";
const string TEST_DATA_PATH="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/dataset_02/";
const string TEST_LIG_DOF_FILE_NAME=TEST_DATA_PATH+"dofLigand.dat";
const string TEST_REC_DOF_FILE_NAME=TEST_DATA_PATH+"dofReceptor.dat";
const string TEST_LIG_POSITION_FILE_NAME=TEST_DATA_PATH+"deformBeforeLigand.dat";
const string TEST_LIG_TRANSFORMED_POSITION_FILE_NAME=TEST_DATA_PATH+"deformAfterLigand.dat";
const string TEST_LIG_MODE_FILE_NAME=TEST_DATA_PATH+"writtenLigandModes.dat";
const string TEST_LIG_FORCE_FILE_NAME=TEST_DATA_PATH+"forcesLigand.dat";
const string TEST_CONFIG_FILE_NAME=TEST_DATA_PATH+"config.dat";

#endif /* TEST_REFERENCEDATA_H_ */
