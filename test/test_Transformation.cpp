#include "test_Transformation.h"
#include <gmock/gmock.h>
#include "test_HMMread.h"
#include "transform.h"
#include "test_ReferenceData.h"
#include "CompareData.h"
#include "Vec3.h"
#include "DOFTransform.h"
#include "Protein.h"
#include "readFile.h"

namespace as {

// this test uses a dataset (DOFS and coordinates before and after deformation and transformation)
//from the CPU version of attract.
/**
 * it then compares the result of the new rotate_translate function to the reference values of the old version
 */
TEST(Transformation_check,compares_if_data_is_transformed_and_deformed_right) {
	//preparation of data
	int numModes=5;
	using namespace std;
	DOF_6D_Modes<float> DofLig;
	DOF_6D_Modes<float> DofRec;
	testProteinConfig<float> testConfig;

	bool centeredReceptor= false;
	bool centeredLigand= true;


	float* inputposX=new float[testConfig.numAtomsLigand];
	float* inputposY=new float[testConfig.numAtomsLigand];
	float* inputposZ=new float[testConfig.numAtomsLigand];

	test_readConfig<float>( &testConfig,TEST_CONFIG_FILE_NAME);
	//std::shared_ptr<Protein<float>> testLigand = createProteinFromPDB<float>("/home/glenn/Documents/Masterthesis/testfolder/1AVX/input/ligandr.pdb") ;

	float* referenceModes= new float[3*testConfig.numAtomsLigand*numModes];
	test_readReferencePositions<float>(inputposX, inputposY, inputposZ ,TEST_LIG_POSITION_FILE_NAME, testConfig.numAtomsLigand);
	float epsilon=0.01;
	CompareData<float> PosX(testConfig.numAtomsLigand,epsilon);
	CompareData<float> PosY(testConfig.numAtomsLigand,epsilon);
	CompareData<float> PosZ(testConfig.numAtomsLigand,epsilon);

	test_readDOF(&DofLig, TEST_LIG_DOF_FILE_NAME, numModes);
	test_readDOF(&DofRec, TEST_REC_DOF_FILE_NAME, numModes);
	DofLig.pos=DofLig.pos+testConfig.pivotLigand;
	std::vector<DOF_6D_Modes<float>> testDofLig;
	testDofLig.push_back(DofLig);
	std::vector<DOF_6D_Modes<float>> testDofRec;
	testDofRec.push_back(DofRec);



	//testLigand->setNumModes(numModes);
	//float* inputPos=testLigand->getOrCreatePosPtr();

	test_readReferenceModes<float>(referenceModes,  TEST_LIG_CONVERTEDMODE_FILE_NAME,  testConfig.numAtomsLigand, numModes);
	test_readReferencePositions<float>(PosX.referenceData(), PosY.referenceData(), PosZ.referenceData(), TEST_LIG_TRANSFORMED_POSITION_FILE_NAME, testConfig.numAtomsLigand);



	//testing transformation of coordinates
	//testLigand->auto_pivotize();
	//transformDOF_glob2rec(testDofRec, testDofLig, testConfig.pivotReceptor, testConfig.pivotLigand, centeredReceptor, centeredLigand);
	rotate_translate(
			inputposX,
			inputposY,
			inputposZ,
		testDofLig.at(0).pos,
		testDofLig.at(0).ang,
		testConfig.numAtomsLigand,
		numModes,
		testDofLig.at(0).modes,
		referenceModes,
		referenceModes+numModes*testConfig.numAtomsLigand,
		referenceModes+2*numModes*testConfig.numAtomsLigand,
		PosX.testData(),
		PosY.testData(),
		PosZ.testData());

	//evaluation
	PosX.evaluateDifference();
	PosY.evaluateDifference();
	PosZ.evaluateDifference();

	PosX.writeResultToFile(EVALUATION_DATA_PATH+"transform_posx.dat");
	PosY.writeResultToFile(EVALUATION_DATA_PATH+"transform_posy.dat");
	PosZ.writeResultToFile(EVALUATION_DATA_PATH+"transform_posz.dat");

	delete [] referenceModes;

	delete [] inputposX;
	delete [] inputposY;
	delete [] inputposZ;

}







}//end as namespace
