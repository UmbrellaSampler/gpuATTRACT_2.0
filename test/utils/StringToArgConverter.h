/*
 * StringToArgConverter.h
 *
 *  Created on: Jan 21, 2018
 *      Author: uwe
 */

#ifndef TEST_UTILS_STRINGTOARGCONVERTER_H_
#define TEST_UTILS_STRINGTOARGCONVERTER_H_

#include <string>
#include <memory>
#include <vector>

namespace as {

class StringToArgConverter {
public:
	explicit StringToArgConverter(std::string cmd) {
		this->_cmd = cmd;
		_argc = 0;
		convert();
	}

	~StringToArgConverter() {
		for (char* arg : _args) {
			delete[] arg;
		}
	}

	int argc() {
		return _argc;
	}

	char** argv() {
		return _args.data();
	}

private:

	void convert() {
		std::istringstream iss(_cmd);
		std::string token;
		while (iss >> token) {
			char *arg = new char[token.size() + 1];
			std::copy(token.begin(), token.end(), arg);
			arg[token.size()] = '\0';
			_args.push_back(arg);
		}
		_args.push_back(nullptr);
		_argc = _args.size() - 1;
	}

	std::string _cmd;
	int _argc;
	std::vector<char *> _args;
};

}



#endif /* TEST_UTILS_STRINGTOARGCONVERTER_H_ */
