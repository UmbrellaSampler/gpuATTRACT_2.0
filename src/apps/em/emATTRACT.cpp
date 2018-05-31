/*
 * emATTRACT.cpp
 *
 *  Created on: Jul 18, 2017
 *      Author: uwe
 */

#include "emATTRACT.tpp"
#include "Chunk.h"
#include "Types_6D.h"

namespace as {

template
class emATTRACT<Types_6D<float>>;

template
class emATTRACT<Types_6D<double>>;

template
class emATTRACT<Types_6D_Modes<float>>;

template
class emATTRACT<Types_6D_Modes<double>>;

}  // namespace as


