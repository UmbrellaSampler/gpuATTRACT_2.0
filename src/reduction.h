/*
 * reduction.h
 *
 *  Created on: Aug 13, 2016
 *      Author: uwe
 */

#ifndef SRC_REDUCTION_H_
#define SRC_REDUCTION_H_

#include "Vec3.h"
#include "Torque.h"
#include "TorqueMat.h"
#include "MatrixFunctions.h"
#include "nativeTypesWrapper.h"
#include "Types_6D.h"

namespace as {

template<typename REAL>
class PotForce {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
public:
	real_t E;
	vec3_t pos;
};


constexpr float ForceLim = 1.0e18;

template<typename REAL>
PotForce<REAL> reducePotForce(
		REAL const* fx, REAL const* fy, REAL const* fz,
		REAL const* energy,
		unsigned const& numAtoms)
{
	PotForce<REAL> potForce;
	potForce.E = 0.0;
	potForce.pos = Vec3<REAL>(0.0);

	typename Types_6D<REAL>::Result res;
	res.E = potForce.E;
	res.pos = potForce.pos;

	for(unsigned i = 0; i < numAtoms; ++i) {
		potForce.E += energy[i];
		potForce.pos.x += fx[i];
		potForce.pos.y += fy[i];
		potForce.pos.z += fz[i];
	}

	// force reduction, some times helps in case of very "bad" start structure
	// taken from original ATTRACT code in trans.f
	for(unsigned i = 0; i < 3; ++i) {
		REAL magn2 = potForce.pos.x*potForce.pos.x
				+ potForce.pos.y*potForce.pos.y
				+ potForce.pos.z*potForce.pos.z;

		if(magn2 > static_cast<REAL>(ForceLim)) {
			potForce.pos.x *= 0.01;
			potForce.pos.y *= 0.01;
			potForce.pos.z *= 0.01;
		}
	}

	return potForce;
}

template<typename REAL>
Vec3<REAL> reduceTorque(
		REAL const* x, REAL const* y, REAL const* z, // unrotated protein coordinates x,y,z !!!
		REAL const* fx, REAL const* fy, REAL const* fz,
		unsigned const& numAtoms,
		Vec3<REAL> const& ang)
{
	Torque<REAL> torque(0);
	for (unsigned i = 0; i < numAtoms; ++i) {
		torque.mat[0][0] += x[i]*fx[i];
		torque.mat[0][1] += y[i]*fx[i];
		torque.mat[0][2] += z[i]*fx[i];
		torque.mat[1][0] += x[i]*fy[i];
		torque.mat[1][1] += y[i]*fy[i];
		torque.mat[1][2] += z[i]*fy[i];
		torque.mat[2][0] += x[i]*fz[i];
		torque.mat[2][1] += y[i]*fz[i];
		torque.mat[2][2] += z[i]*fz[i];
	}

	const TorqueMat<REAL> torqueMat = euler2torquemat(ang.x, ang.y, ang.z);
	Vec3<REAL> result = torqueMat.rotateReduce(torque);

	return result;
}

inline bool isPow2(unsigned int x)
{
    return ((x&(x-1))==0);
}

inline unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

#ifdef CUDA

template <class T>
void d_reduce(
		const unsigned& threads,
		const unsigned& blocks,
		const unsigned& numAtoms,
		T* xPos, T* yPos, T* zPos,
		T *d_fx, T *d_fy, T *d_fz, T *d_E,
		T *g_odata,
		const cudaStream_t& stream);

/* wrapper function to call the device kernel for the reduction */
template<typename REAL>
void deviceReduce(
		const unsigned& blockSize,
		const unsigned& numDOFs,
		const unsigned& ligandSize,
		REAL* xPos, REAL* yPos, REAL* zPos,
		REAL *d_fx, REAL *d_fy, REAL *d_fz, REAL *d_E,
		REAL *d_out,
		const cudaStream_t& stream)
{
	/* we need at least twice the block size number of threads */
	const unsigned threads = (ligandSize < blockSize*2) ? nextPow2((ligandSize + 1)/ 2) : blockSize;
	/* each structure is reduced by one thread block */
	const unsigned blocks = numDOFs;
	d_reduce(threads, blocks, ligandSize,
			xPos, yPos, zPos,
			d_fx, d_fy, d_fz, d_E,
			d_out,
			stream);
}

/* remaining reduce part that is performed on the host */
template<typename REAL>
void h_finalReduce(
			const unsigned& numDOFs,
			typename Types_6D<REAL>::DOF* dofs,
			const REAL* deviceOut,
			typename Types_6D<REAL>::Result* enGrads)
{

	for (unsigned i = 0; i < numDOFs; ++i)
	{
		auto &enGrad = enGrads[i];
		enGrad.pos.x = deviceOut[i*13 + 0];
		enGrad.pos.y = deviceOut[i*13 + 1];
		enGrad.pos.z = deviceOut[i*13 + 2];

		for(unsigned j = 0; j < 3; ++j) {
			REAL magn2 = enGrad.pos.x*enGrad.pos.x
					+ enGrad.pos.y*enGrad.pos.y
					+ enGrad.pos.z*enGrad.pos.z;

			if(magn2 > static_cast<REAL>(ForceLim)) {
				enGrad.pos.x *= 0.01;
				enGrad.pos.y *= 0.01;
				enGrad.pos.z *= 0.01;
			}
		}

		enGrad.E = deviceOut[i*13 + 3];

		Torque<REAL> torque;
		torque.mat[0][0] = deviceOut[i*13 + 4 ];
		torque.mat[0][1] = deviceOut[i*13 + 5 ];
		torque.mat[0][2] = deviceOut[i*13 + 6 ];
		torque.mat[1][0] = deviceOut[i*13 + 7 ];
		torque.mat[1][1] = deviceOut[i*13 + 8 ];
		torque.mat[1][2] = deviceOut[i*13 + 9 ];
		torque.mat[2][0] = deviceOut[i*13 + 10];
		torque.mat[2][1] = deviceOut[i*13 + 11];
		torque.mat[2][2] = deviceOut[i*13 + 12];

		const auto &dof = dofs[i];
		const TorqueMat<REAL> torqueMat = euler2torquemat(dof.ang.x, dof.ang.y, dof.ang.z);
		Vec3<REAL> result = torqueMat.rotateReduce(torque);

		enGrad.ang.x = result.x;
		enGrad.ang.y = result.y;
		enGrad.ang.z = result.z;
	}
}

#endif // CUDA


} // namespace as


#endif /* SRC_REDUCTION_H_ */
