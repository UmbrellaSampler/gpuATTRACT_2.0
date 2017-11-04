#include "Types_6D_modes.tpp"

using namespace std;

namespace as {

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D_Modes<float> const& dof);

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D_Modes<double> const& dof);


template
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<float> const& enGrad);

template
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<double> const& enGrad);

}  // namespace as

