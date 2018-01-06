#include "Types_6D_MB_Modes.tpp"

using namespace std;

namespace as {

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D_MB_Modes<float> const& dof);

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D_MB_Modes<double> const& dof);


template
std::ostream& operator<< (std::ostream& s, Result_6D_MB_Modes<float> const& enGrad);

template
std::ostream& operator<< (std::ostream& s, Result_6D_MB_Modes<double> const& enGrad);

unsigned int Common_MB_Modes::numModesRec = 0;
unsigned int Common_MB_Modes::numModesLig[LIGANDS_MAX_NUMBER] = {0};
unsigned int Common_MB_Modes::numLigands = 1;
}  // namespace as

