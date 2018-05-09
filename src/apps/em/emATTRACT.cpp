/*
 * emATTRACT.cpp
 *
 *  Created on: Jul 18, 2017
 *      Author: uwe
 */

#include "emATTRACT.tpp"
#include "Chunk.h"
#include "Types_6D.h"
#include "Types_6D_Modes.h"

namespace as {

template
class emATTRACT<Types_6D<float>>;

template
class emATTRACT<Types_6D<double>>;

// TODO: is implemented after scATTRACT<Types_6D_Modes<>>
template
class emATTRACT<Types_6D_Modes<float>>;

template
class emATTRACT<Types_6D_Modes<double>>;



}  // namespace as


