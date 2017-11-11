/*
 * test_HMMread.h
 *
 *  Created on: Nov 11, 2017
 *      Author: glenn
 */

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


#ifndef TEST_HMMREAD_H_
#define TEST_HMMREAD_H_

namespace as {


template<typename REAL>
struct testProteinConfig{
	int numAtomsReceptor;
	int numAtomsLigand;
	Vec3<REAL> pivotReceptor;
	Vec3<REAL> pivotLigand;
};

static std::vector<std::string> line2Strings(std::string const& line) {
	using namespace std;
	istringstream iss(line);
	return vector<string> { istream_iterator<string> { iss },
					istream_iterator<string> { } };
}


template<typename REAL>
void readConfig( testProteinConfig<REAL> *testconfig,std::string configFileName){
	using namespace std;
	fstream configFile(configFileName);
	string line;
	vector<std::string> tokens;
	bool receptor=false;
	bool natom=false;

	if (!configFile.is_open()){	perror(("error while opening file " + configFileName).c_str());}
	while(getline(configFile, line)) {
		if(!line.empty()){



			tokens=line2Strings(line);

			if(tokens.at(0).compare("##")==0){

				if (tokens.at(1).compare("receptor")==0){receptor=true;}
				else if( tokens.at(1).compare("ligand")==0){receptor=false;}
				else if (tokens.at(1).compare("natom")==0){natom=true;}
				else if (tokens.at(1).compare("pivot")==0){natom=false;}
				}
			if(tokens.at(0).compare("#")==0){

				if (receptor && natom){testconfig->numAtomsReceptor=stoi(tokens.at(1));}
				else if (receptor && !natom){testconfig->pivotReceptor.x=stof(tokens.at(1));testconfig->pivotReceptor.y=stof(tokens.at(2));testconfig->pivotReceptor.z=stof(tokens.at(3));}
				else if (!receptor && natom){testconfig->numAtomsLigand=stoi(tokens.at(1));}
				else if (!receptor && !natom){testconfig->pivotLigand.x=stof(tokens.at(1));testconfig->pivotLigand.y=stof(tokens.at(2));testconfig->pivotLigand.z=stof(tokens.at(3));}
			}
		}
	}
	std::cout<<testconfig->numAtomsReceptor<<"\n"<<testconfig->numAtomsLigand<<"\n"<<testconfig->pivotLigand.x<<"\n"<<testconfig->pivotLigand.y<<"\n"<<testconfig->pivotLigand.z<<std::endl;
}






template<typename REAL>
void test_readReferenceModes(REAL *referenceModes, std::string modeFileName, int numAtoms, int numModes) {
	using namespace std;

	vector<REAL> modeX;	vector<REAL> modeY;	vector<REAL> modeZ;
	REAL modeVal[numModes];
	REAL eigVal[numModes];
	protModes = new REAL[3*numModes*numAtoms];

	bool changeMode=true;
	bool isData=false;
	int idxPos=0;
	int idxLine=0;
	int modeIdx=0;

	vector<std::string> tokens;
	string line;
	ifstream modeFile(modeFileName);

	if (!modeFile.is_open()){	perror(("error while opening file " + modeFileName).c_str());}

	while(getline(modeFile, line)) {
		if(!line.empty()){
			tokens=line2Strings(line);

			if(idxPos < numAtoms && changeMode==false){isData=true;}
			else if(idxPos == numAtoms){isData=false;idxPos=0;}

			if(changeMode == true && tokens.size() > 1){
				modeIdx++;
				if (modeIdx == stoi(tokens.at(0)) && (modeIdx <= numModes)) {
					modeVal[modeIdx]=stof(tokens.at(1))*stof(tokens.at(1));
					eigVal[modeIdx-1]=0;
					//cout <<"# "<<modeVal[modeIdx]<< endl;
					changeMode=false;
				}
			}
			if(isData==true && tokens.size() == 3 ){
				float x=stof(tokens.at(0));
				float y=stof(tokens.at(1));
				float z=stof(tokens.at(2));

				eigVal[modeIdx-1]+=x*x+y*y+z*z;
				referenceModes[numModes*idxPos+modeIdx-1]						=stof(tokens.at(0));
				referenceModes[numAtoms*numModes+numModes*idxPos+modeIdx-1]		=stof(tokens.at(1));
				referenceModes[2*numAtoms*numModes+numModes*idxPos+modeIdx-1]	=stof(tokens.at(2));
				idxPos++;
				if(idxPos==numAtoms){changeMode=true;}
			}
		}
		idxLine++;
	}
	if (modeFile.bad()){	perror(("error while reading file " + modeFileName).c_str());}
	modeFile.close();

}

template
void test_readReferenceModes(float *referenceModes, std::string modeFileName, int numAtoms, int numModes);

template
void test_readReferenceModes(double *referenceModes, std::string modeFileName, int numAtoms, int numModes);

template
void readConfig( testProteinConfig<float> *testconfig,std::string configFileName);

template
void readConfig( testProteinConfig<double> *testconfig,std::string configFileName);

template
struct testProteinConfig<float>;

template
struct testProteinConfig<double>;
}
#endif /* TEST_HMMREAD_H_ */
