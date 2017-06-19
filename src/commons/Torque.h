/*
 * Torque.h
 *
 *  Created on: Aug 13, 2016
 *      Author: uwe
 */

#ifndef SRC_TORQUE_H_
#define SRC_TORQUE_H_

namespace as {

template<typename T>
class Torque {
public:
	Torque() {};
	Torque(T value): mat{value} {}

	T mat[3][3];
};

}  // namespace as



#endif /* SRC_TORQUE_H_ */
