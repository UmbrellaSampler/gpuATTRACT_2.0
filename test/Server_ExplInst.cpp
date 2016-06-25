/*
 * Server_FwdDecl.cpp
 *
 *  Created on: Feb 6, 2016
 *      Author: uwe
 */

#include "../src/ServerIncludes.h"

#include "Service_Mock.h"
#include "Service_TimeOut.h"

template
class as::Server<test::Service_Mock>;

template
class as::Server<test::Service_TimeOut>;

template
class as::HostAllocator<test::Service_Mock::input_t>;
