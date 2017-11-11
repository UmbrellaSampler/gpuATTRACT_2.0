#include "test_Transformation.h"
#include <gmock/gmock.h>
#include "test_HMMread.h"
#include "transform.h"

namespace as {

TEST(a,s) {


DOF_6D_Modes<float> testDof;
std::string dofFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/1_dofLigand.dat";
string configFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/1_config.dat";
string referenceModeFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/reference/1_writtenLigandModes.dat";
string modeFileName="/home/glenn/Documents/Masterthesis/testfolder/1AVX/input/modesLigand.dat";

test_readDOF<float>( &testDof,dofFileName, 5);

std::cout<<testDof<<std::endl;




//	rotate_translate<float>(
//			REAL const* x,
//			REAL const* y,
//			REAL const* z,
//			Vec3<REAL> const& displacement,
//			Vec3<REAL> const& ang,
//			unsigned const& ligSize,
//			unsigned const& numModes,
//			REAL const* dlig,
//			REAL const* xModes,
//			REAL const* yModes,
//			REAL const* zModes,
//			REAL* xTr,
//			REAL* yTr,
//			REAL* zTr)

}






}//end as namespace
