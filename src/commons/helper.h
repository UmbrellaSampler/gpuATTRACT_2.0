/*
 * helper.h
 *
 *  Created on: Mar 21, 2016
 *      Author: uwe
 */

#ifndef SRC_HELPER_H_
#define SRC_HELPER_H_

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


#endif /* SRC_HELPER_H_ */
