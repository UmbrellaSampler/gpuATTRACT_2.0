/*
 * test_Transformation.h
 *
 *  Created on: Nov 11, 2017
 *      Author: glenn
 */

#ifndef TEST_TRANSFORMATION_H_
#define TEST_TRANSFORMATION_H_

#include "readFile.h"
#include <sstream>
#include <fstream>
#include <cstring>
#include "Protein.h"
#include <iostream>
#include <iterator>
#include "Vec3.h"
#include <string>
#include <memory>
#include "Types_6D_Modes.h"



static std::vector<std::string> line2Strings(std::string const& line) {
	using namespace std;
	istringstream iss(line);
	return vector<string> { istream_iterator<string> { iss },
					istream_iterator<string> { } };
}

namespace as {



//template<typename REAL>
//void test_readReferencePositions(REAL* referencePosX,REAL* referencePosY,REAL* referencePosZ, std::string positionFileName, int numAtoms) {
//	using namespace std;
//
//	vector<REAL> posX;	vector<REAL> posY;	vector<REAL> posZ;
//
//	vector<std::string> tokens;
//	string line;
//	ifstream positionFile(positionFileName);
//	int idx=0;
//	if (!positionFile.is_open()){	perror(("error while opening file " + positionFileName).c_str());}
//
//	while(getline(positionFile, line)) {
//		if(!line.empty()){
//			tokens=line2Strings(line);
//			referencePosX[idx]=stof(tokens.at(0));
//			referencePosY[idx]=(stof(tokens.at(1)));
//			referencePosZ[idx]=(stof(tokens.at(2)));
//			idx++;
//			}
//		}
//	if (positionFile.bad()){	perror(("error while reading file " + positionFileName).c_str());}
//	positionFile.close();
//}
//
//
//template<typename REAL>
//void test_readDOF( DOF_6D_Modes<REAL> *testDof,std::string dofFileName, int numModes){
//	using namespace std;
//	fstream dofFile(dofFileName);
//	string line;
//	vector<std::string> tokens;
//
//	bool read=false;
//	int dofidx=0;
//	testDof->numModes=5;
//	if (!dofFile.is_open()){	perror(("error while opening file " + dofFileName).c_str());}
//	while(getline(dofFile, line)) {
//		if(!line.empty()){
//
//			tokens=line2Strings(line);
//			if(tokens.at(0).compare("#")==0){dofidx++;read=true;}
//			if(tokens.at(0).compare("###")==0){read=false;}
//			if (dofidx==1 && read==true){for(int mode=0;mode<numModes;mode++){testDof->modes[mode]=stof(tokens.at(mode+1));}}
//			if (dofidx==2 && read==true){testDof->pos.x=stof(tokens.at(1));testDof->pos.y=stof(tokens.at(2));testDof->pos.z=stof(tokens.at(3));}
//			if (dofidx==3 && read==true){testDof->ang.x=stof(tokens.at(1));testDof->ang.y=stof(tokens.at(2));testDof->ang.z=stof(tokens.at(3));}
//
//		}
//		}
//	if (dofFile.bad()){	perror(("error while reading file " + dofFileName).c_str());}
//	dofFile.close();
//}


//template
//void test_readReferencePositions(float* referencePosX,float* referencePosY,float* referencePosZ, std::string positionFileName, int numAtoms) ;
//
//
//template
//void test_readReferencePositions(double* referencePosX,double* referencePosY,double* referencePosZ, std::string positionFileName, int numAtoms) ;
//
//
//
//template
//void test_readDOF( DOF_6D_Modes<float> *testDof,std::string dofFileName, int numModes);
//
//template
//void test_readDOF( DOF_6D_Modes<double> *testDof,std::string dofFileName, int numModes);

}
#endif /* TEST_TRANSFORMATION_H_ */
