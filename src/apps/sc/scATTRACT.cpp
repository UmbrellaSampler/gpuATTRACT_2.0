/*
 * scATTRACT.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#include "scATTRACT.tpp"

namespace as {

template
class scATTRACT<Types_6D<float>>;

template
class scATTRACT<Types_6D<double>>;

template
class scATTRACT<Types_6D_Modes<float>>;

template
class scATTRACT<Types_6D_Modes<double>>;

template
class scATTRACT<Types_6D_MB_Modes<float>>;

template
class scATTRACT<Types_6D_MB_Modes<double>>;

}  // namespace as

