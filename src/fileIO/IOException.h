/*
 * IOException.h
 *
 *  Created on: Jan 13, 2018
 *      Author: uwe
 */

#ifndef SRC_FILEIO_IOEXCEPTION_H_
#define SRC_FILEIO_IOEXCEPTION_H_

#include <stdexcept>
#include <string>

namespace as {

class IOException : public std::logic_error {
public:
	explicit IOException (const std::string& what_arg) :
			logic_error(what_arg) {}
	explicit IOException (const char* what_arg) :
			logic_error(what_arg) {}
};

}  // namespace as



#endif /* SRC_FILEIO_IOEXCEPTION_H_ */
