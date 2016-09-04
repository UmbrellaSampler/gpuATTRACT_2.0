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

#ifndef UTILS_H_
#define UTILS_H_

#include <cstdio>
#include <cstdlib>

#ifdef CUDA
#include <cuda_runtime.h>

#ifndef NDEBUG
#define cudaVerify(x) do { 																				\
		cudaError_t __cu_result = x; 																	\
		if (__cu_result!=cudaSuccess) { 																\
			fprintf(stderr, "%s:%i: Error: cuda function call failed:\n" 								\
					"%s;\nmessage: %s\n", 																\
					__FILE__, __LINE__, #x, cudaGetErrorString(__cu_result));							\
			exit(1);																					\
		} 																								\
	} while(0)
#define cudaVerifyKernel(x) do {																		\
		x;																								\
		cudaDeviceSynchronize();																		\
		cudaError_t __cu_result = cudaGetLastError();													\
		if (__cu_result!=cudaSuccess) { 																\
			fprintf(stderr, "%s:%i: Error: cuda function call failed:\n" 								\
					"%s;\nmessage: %s\n", 																\
					__FILE__, __LINE__, #x, cudaGetErrorString(__cu_result));							\
			exit(1);																					\
		} 																								\
	} while(0)
#else
#define cudaVerify(x) do {																				\
		x;																								\
	} while(0)
#define cudaVerifyKernel(x) do {																		\
		x;																								\
	} while(0)
#endif

#define CUDA_CHECK(x) do {                                                        \
	cudaError_t __cu_result = x;                                                  \
	if (__cu_result!=cudaSuccess) { 										      \
		char buffer[1000];                                                        \
		sprintf(buffer, "%s:%i: Error: cuda function call failed:\n" 		      \
				"%s;\nmessage: %s\n", 										      \
				__FILE__, __LINE__, #x, cudaGetErrorString(__cu_result));         \
		throw std::runtime_error(std::string(buffer));                            \
		exit(1);															      \
	} 																		      \
} while(0)

#endif

#define ASSERT(x) do {                                           \
	if((x) == false) {                                           \
		fprintf(stderr, "%s:%i: Assertion '%s' failed.\n",     	\
				__FILE__, __LINE__, #x );                        \
		exit(1);                                                 \
	}                                                            \
} while(0)




#endif /* UTILS_H_ */
