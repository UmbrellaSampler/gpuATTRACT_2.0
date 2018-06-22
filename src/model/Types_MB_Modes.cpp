#include "Types_MB_Modes.tpp"

using namespace std;

namespace as {

template
std::ostream& operator <<(std::ostream& outStream, DOF_MB_Modes<float> const& dof);

template
std::ostream& operator <<(std::ostream& outStream, DOF_MB_Modes<double> const& dof);


template
std::ostream& operator<< (std::ostream& s, Result_MB_Modes<float> const& enGrad);

template
std::ostream& operator<< (std::ostream& s, Result_MB_Modes<double> const& enGrad);

unsigned int Common_MB_Modes::numModes[NUM_MAX_PROTEIN]={0};


}  // namespace as

