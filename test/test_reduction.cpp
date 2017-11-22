#include "test_Transformation.h"
#include "test_reduction.h"

namespace as{


TEST(Reduction,checks_for_right_reduction_of_allforces){
	int numModes=5;
	testProteinConfig<float> testConfig;
	test_readConfig<float>( &testConfig,TEST_CONFIG_FILE_NAME);
	DOF_6D_Modes<float> DofLig;

	float* forceX=new float[testConfig.numAtomsLigand];
	float* forceY=new float[testConfig.numAtomsLigand];
	float* forceZ=new float[testConfig.numAtomsLigand];

	float* posX=new float[testConfig.numAtomsLigand];
	float* posY=new float[testConfig.numAtomsLigand];
	float* posZ=new float[testConfig.numAtomsLigand];

	float* origposX=new float[testConfig.numAtomsLigand];
	float* origposY=new float[testConfig.numAtomsLigand];
	float* origposZ=new float[testConfig.numAtomsLigand];

	float* modes=new float[3*numModes*testConfig.numAtomsLigand];


	test_readReferencePositions<float>(posX, posY, posZ, TEST_LIG_TRANSFORMED_POSITION_FILE_NAME, testConfig.numAtomsLigand);
	test_readReferencePositions<float>(origposX, origposY, origposZ, TEST_LIG_POSITION_FILE_NAME, testConfig.numAtomsLigand);
	test_readReferenceModes<float>(modes,  TEST_LIG_CONVERTEDMODE_FILE_NAME,  testConfig.numAtomsLigand, numModes);

	test_readDOF(&DofLig, TEST_LIG_DOF_FILE_NAME, numModes);
	test_readReferenceForce(forceX,forceY,forceZ, TEST_LIG_FORCE_FILE_NAME, testConfig.numAtomsLigand);

	PotForce<float> redPotForce = reducePotForce(
			forceX,
			forceY,
			forceZ,
			forceX,
			testConfig.numAtomsLigand
	);


	reduceModeForce(
			DofLig.ang,
			forceX,
			forceY,
			forceZ,
			modes,
			modes+numModes*testConfig.numAtomsLigand,
			modes+numModes*testConfig.numAtomsLigand*2,
			testConfig.numAtomsLigand,
			numModes,
			redPotForce.modes
			);

	Result_6D_Modes<float> enGrad;
	enGrad.numModes=numModes;
	for( int mode=0;mode<numModes;mode++){

		enGrad.modes[mode]=redPotForce.modes[mode];}
	enGrad.E = redPotForce.E;
	enGrad.pos = redPotForce.pos;

	enGrad.ang = reduceTorque(
			origposX,
			origposY,
			origposZ,
			forceX,
			forceY,
			forceZ,
			testConfig.numAtomsLigand,
			DofLig.ang
	); // OK

	//print Gradients
	std::cout<<enGrad<<std::endl;


	//detete data
	delete [] forceX;
	delete [] forceY;
	delete [] forceZ;

	delete [] posX;
	delete [] posY;
	delete [] posZ;

	delete [] origposX;
	delete [] origposY;
	delete [] origposZ;

	delete [] modes;
}




}
