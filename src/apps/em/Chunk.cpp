/*
 * Chunk.cpp
 *
 *  Created on: Jul 9, 2017
 *      Author: uwe
 */

#include "Chunk.tpp"
#include "Server.h"
#include "Request.h"

namespace as {

template
class RequestHandler<Types_6D<float>>::Chunk;

template
class RequestHandler<Types_6D<double>>::Chunk;

template
class RequestHandler<Types_6D_Modes<float>>::Chunk;

template
class RequestHandler<Types_6D_Modes<double>>::Chunk;

template
class RequestHandler<Types_MB_Modes<float>>::Chunk;

template
class RequestHandler<Types_MB_Modes<double>>::Chunk;

}  // namespace as


