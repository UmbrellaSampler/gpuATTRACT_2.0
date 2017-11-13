#include <gmock/gmock.h>
#include "test_HMMread.h"

namespace as {


TEST(ReadModefunction, ckecksmodeconversion) {

	using namespace std;
	testProteinConfig<float> testConfig;
	string configFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/1_config.dat";
	string referenceModeFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/1_writtenLigandModes.dat";
	string modeFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/input/modesLigand.dat";

	test_readConfig( &testConfig, configFileName);

	float* referenceModes;
	int numModes=5;
	float e=0.01;
	float xmode,ymode,zmode;
	int count=0,c1=0;
	referenceModes=(float*) malloc(3*testConfig.numAtomsLigand*numModes*sizeof(float));
	test_readReferenceModes(referenceModes,  referenceModeFileName, testConfig.numAtomsLigand, numModes);
	//std::cout<<"test"<<std::endl;

	std::shared_ptr<Protein<float>> testProt = std::make_shared<Protein<float>>();
	testProt->setNumAtoms(testConfig.numAtomsLigand);
	testProt->setNumModes(numModes);

	readHMMode(testProt,modeFileName);


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

	//std::cout<< 3*c1<<" "<<count<<std::endl;
	free(referenceModes);
}


}// end namespace as

