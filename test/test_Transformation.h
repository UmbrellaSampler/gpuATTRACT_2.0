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



template<typename REAL>
void test_readReferencePosiotions(REAL* referencePosX,REAL* referencePosY,REAL* referencePosZ, std::string positionFileName, int numAtoms) {
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
			referencePosX[idx]=stof(tokens.at(0));
			referencePosY[idx]=(stof(tokens.at(1)));
			referencePosZ[idx]=(stof(tokens.at(2)));
			idx++;
			}
		}
	if (positionFile.bad()){	perror(("error while reading file " + positionFileName).c_str());}
	positionFile.close();
}


template<typename REAL>
void test_readDOF( DOF_6D_Modes<REAL> testDof,std::string dofFileName, int numModes){
	using namespace std;
	fstream dofFile(dofFileName);
	string line;
	vector<std::string> tokens;
	bool receptor=false;
	bool natom=false;
	int dofidx=0;
	if (!dofFile.is_open()){	perror(("error while opening file " + dofFileName).c_str());}
	while(getline(dofFile, line)) {
		if(!line.empty()){

			tokens=line2Strings(line);
			if(tokens.at(0).compare("#")==0){dofidx++;}
			if (dofidx==1){for(int mode=0;mode<numModes;mode++){testDof->modes[mode]=tokens.at(mode+1);}}
			if (dofidx==2){testDof->pos.x=tokens.at(1);testDof->pos.y=tokens.at(2);testDof->pos.z=tokens.at(3);}
			if (dofidx==3){testDof->ang.x=tokens.at(1);testDof->ang.y=tokens.at(2);testDof->ang.z=tokens.at(3);}
			}
		}
	if (dofFile.bad()){	perror(("error while reading file " + dofFileName).c_str());}
	dofFile.close();
}

}
#endif /* TEST_TRANSFORMATION_H_ */