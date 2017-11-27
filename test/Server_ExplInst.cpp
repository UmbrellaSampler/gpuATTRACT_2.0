/*
 * Server_FwdDecl.cpp
 *
 *  Created on: Feb 6, 2016
 *      Author: uwe
 */

#include "ServerIncludes.h"

#include "Service_Mock.h"
#include "Service_TimeOut.h"

template
class as::Server<test::Service_Mock::genericTypes_t>;

//template
//class as::Server<test::Service_TimeOut::genericTypes_t>;

template
class as::HostAllocator<test::Service_Mock::input_t>;
