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


namespace as {

template
class RequestHandler<Types_6D<float>>;

template
class RequestHandler<Types_6D<double>>;

template
class RequestHandler<Types_6D_Modes<float>>;

template
class RequestHandler<Types_6D_Modes<double>>;

}  // namespace as

