/*
 * ServerBuffers.h
 *
 *  Created on: Jan 1, 2016
 *      Author: uwe
 */

#ifndef SERVERBUFFERS_H_
#define SERVERBUFFERS_H_

namespace as {

template<typename InputType, typename ResultType>
struct ServerBuffers {
	InputType* inputBuffer;
	ResultType* resultBuffer;
	size_t size;
};

}  // namespace as



#endif /* SERVERBUFFERS_H_ */
