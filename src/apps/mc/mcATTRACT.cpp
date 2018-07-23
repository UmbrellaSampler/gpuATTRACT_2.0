/*
 * mcATTRACT.cpp
 *
 *  Created on: Aug 16, 2016
 *      Author: uwe
 */

#include "mcATTRACT.tpp"
#include "Chunk.h"
#include "Types_6D.h"

namespace as {

template
class mcATTRACT<Types_6D<float>>;

template
class mcATTRACT<Types_6D<double>>;

template
class mcATTRACT<Types_6D_Modes<float>>;

template
class mcATTRACT<Types_6D_Modes<double>>;

}  // namespace as

