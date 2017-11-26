#include <gmock/gmock.h>
#include "test_HMMread.h"
#include "test_ReferenceData.h"
#include "CompareData.h"

namespace as {


// this TEST checks weather the modes of a Protein are read and converted just like in the old attract version
TEST(ReadModefunction, ckecksmodeconversion) {
	//data is read and prepared such it can be fed to readHMMmode()

	using namespace std;
	int numModes=5;
	float e=0.01;
	float* referenceModes;
	float xmode,ymode,zmode;
	int count=0,c1=0;

	testProteinConfig<float> testConfig;

	string configFileName=TEST_CONFIG_FILE_NAME;
	string referenceModeFileName=TEST_LIG_CONVERTEDMODE_FILE_NAME;
	string modeFileName=TEST_LIG_MODE_FILE_NAME;

	test_readConfig( &testConfig, configFileName);

	referenceModes=(float*) malloc(3*testConfig.numAtomsLigand*numModes*sizeof(float));

	test_readReferenceModes(referenceModes,  referenceModeFileName, testConfig.numAtomsLigand, numModes);
	std::shared_ptr<Protein<float>> testProt = std::make_shared<Protein<float>>();
	testProt->setNumAtoms(testConfig.numAtomsLigand);
	testProt->setNumModes(numModes);


	//readHMMode function is tested
	readHMMode(testProt,modeFileName);


	//data is evaluated and printed if not within a cerrain ratio
	int numAtoms=testConfig.numAtomsLigand;
	for(int i=0;i<numAtoms;i++){
		for(int mode=0;mode<numModes;mode++){
			xmode=testProt->xModes()[numModes*i+mode];
			ymode=testProt->yModes()[numModes*i+mode];
			zmode=testProt->zModes()[numModes*i+mode];
			c1++;
			if(1-e < referenceModes[numModes*i+mode]/xmode &&  referenceModes[numModes*i+mode]/xmode < 1+e){count++;}
			else{std::cout<<"i: "<< i<<" mode: "<<mode<<"\t ref: "<<referenceModes[numModes*i+mode]<<"\t test: "<<xmode<<std::endl;}
			if(1-e < referenceModes[numAtoms*numModes+numModes*i+mode]/ymode && referenceModes[numAtoms*numModes+numModes*i+mode]/ymode < 1+e){count++;}
			else{std::cout<<"i: "<< i<<" mode: "<<mode<<"\t ref: "<<referenceModes[numAtoms*numModes+numModes*i+mode]<<"\t test: "<<ymode<<std::endl;}
			if(1-e < referenceModes[2*numAtoms*numModes+numModes*i+mode]/zmode && referenceModes[2*numAtoms*numModes+numModes*i+mode]/zmode < 1+e){count++;}
			else{std::cout<<"i: "<< i<<" mode: "<<mode<<"\t ref: "<<referenceModes[2*numAtoms*numModes+numModes*i+mode]<<"\t test: "<<zmode<<std::endl;}
		}
	}

	free(referenceModes);
}


}// end namespace as

