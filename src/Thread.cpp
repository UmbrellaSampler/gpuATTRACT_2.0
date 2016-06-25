/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#include "Thread.h"

/* Constructor */
as::Thread::Thread() {}

/* Destructor */
as::Thread::~Thread() {}


/****************************
 * public member functions
 ****************************/
std::thread::id as::Thread::id() {
	return _thread.get_id();
}
void as::Thread::start() {
	_thread = std::thread(&Thread::run, this);
}
void as::Thread::join() {
	_thread.join();
}
void as::Thread::detach() {
	_thread.detach();
}
void as::Thread::joinable() {
	_thread.joinable();
}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


