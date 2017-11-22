/*
 * test_reduction.h
 *
 *  Created on: Nov 17, 2017
 *      Author: glenn
 */

#ifndef TEST_REDUCTION_H_
#define TEST_REDUCTION_H_

#include <gmock/gmock.h>
#include "test_HMMread.h"
#include "transform.h"
#include "test_ReferenceData.h"
#include "CompareData.h"
#include "Vec3.h"
#include "DOFTransform.h"
#include "Protein.h"
#include "readFile.h"
#include "reduction_modes.h"

template<typename REAL>
void test_readReferenceForce(REAL* referenceForceX,REAL* referenceForceY,REAL* referenceForceZ, std::string positionFileName, int numAtoms) {
	using namespace std;

	vector<REAL> posX;	vector<REAL> posY;	vector<REAL> posZ;

	vector<std::string> tokens;
	string line;
	ifstream positionFile(positionFileName);
	int idx=0;
	if (!positionFile.is_open()){	perror(("error while opening file " + positionFileName).c_str());}

	while(getline(positionFile, line)) {
		if(!line.empty()){
			tokens=line2Strings(line);
			referenceForceX[idx]=stof(tokens.at(0));
			referenceForceY[idx]=(stof(tokens.at(1)));
			referenceForceZ[idx]=(stof(tokens.at(2)));
			idx++;
			}
		}
	if (positionFile.bad()){	perror(("error while reading file " + positionFileName).c_str());}
	positionFile.close();
}

template
void test_readReferenceForce<float>(float* referenceForceX,float* referenceForceY,float* referenceForceZ, std::string positionFileName, int numAtoms);

template
void test_readReferenceForce<double>(double* referenceForceX,double* referenceForceY,double* referenceForceZ, std::string positionFileName, int numAtoms);
#endif /* TEST_REDUCTION_H_ */
