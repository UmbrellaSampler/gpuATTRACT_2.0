/*
 * RequestHandler.cpp
 *
 *  Created on: Jul 9, 2017
 *      Author: uwe
 */


#include "RequestHandler.tpp"
#include "Server.h"
#include "Request.h"
#include "Types_6D.h"
#include "Types_MB_Modes.h"


namespace as {

template
class RequestHandler<Types_6D<float>>;

template
class RequestHandler<Types_6D<double>>;

template
class RequestHandler<Types_6D_Modes<float>>;

template
class RequestHandler<Types_6D_Modes<double>>;

template
class RequestHandler<Types_MB_Modes<float>>;

template
class RequestHandler<Types_MB_Modes<double>>;
}  // namespace as

