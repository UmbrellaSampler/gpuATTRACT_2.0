/*
 * align.h
 *
 *  Created on: Apr 3, 2016
 *      Author: uwe
 */

#ifndef SRC_ALIGN_H_
#define SRC_ALIGN_H_

#if defined(__GNUC__)

#define __align__(n) \
		__attribute__((aligned(n)))

#define __builtin_align__(a) \
        __align__(a)

#elif defined(__CUDACC__)

#include "host_defines.h"

#endif // ifndef CUDA

#endif /* SRC_ALIGN_H_ */
