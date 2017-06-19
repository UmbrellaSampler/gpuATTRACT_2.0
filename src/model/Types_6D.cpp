/*
 * DOF_6D.cpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#include "Types_6D.tpp"

using namespace std;

namespace as {

template
std::ostream& operator <<(std::ostream& outStream, DOFImpl<float> const& dof);

template
std::ostream& operator <<(std::ostream& outStream, DOFImpl<double> const& dof);


template
std::ostream& operator<< (std::ostream& s, ResultImpl<float> const& enGrad);

template
std::ostream& operator<< (std::ostream& s, ResultImpl<double> const& enGrad);

}  // namespace as


