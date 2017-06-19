/*
 * forcefield.h
 *
 *  Created on: Aug 13, 2016
 *      Author: uwe
 */

#ifndef SRC_FORCEFIELD_H_
#define SRC_FORCEFIELD_H_

#include "ParamTable.h"
#include "SimParam.h"

#ifndef __CUDACC__

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

#endif

namespace as {

template<typename REAL>
__host__ __device__
void LJPotForce(
		REAL const& dr2,
		REAL const& dr2_inv,
		REAL const& dx,
		REAL const& dy,
		REAL const& dz,
		typename ParamTable<REAL>::type_t const& LJParams,
		REAL const& swi,
		PotShape const& potShape,
		REAL& fx,
		REAL& fy,
		REAL& fz,
		REAL& energy)
{
	REAL const& alen = LJParams.ac; // 24*eps
	REAL const& rlen = LJParams.rc; // sigma squared
	REAL const& ivor = LJParams.ipon;
	REAL const& rmin2 = LJParams.rmin2;
	REAL const& emin = LJParams.emin;

	REAL const rr23 = dr2_inv*dr2_inv*dr2_inv;

	REAL rrd, shapedelta;

	switch (potShape) {
	case PotShape::_8_6:
		rrd = dr2_inv;
		shapedelta = 2;
		break;
	case PotShape::_12_6:
		rrd = rr23;
		shapedelta = 6;
		break;
	default:
		break;
	}
	REAL rep = rlen * rrd;
	REAL vlj = (rep - alen) * rr23;
	if (dr2 < rmin2) {
		energy = swi * (vlj + (ivor - 1) * emin);
		REAL fb_swi = (6.0 * vlj + shapedelta * (rep * rr23))*swi;
		fx = fb_swi * dx;
		fy = fb_swi * dy;
		fz = fb_swi * dz;

	} else {
		REAL swi_ivor = swi * ivor;
		energy = swi_ivor * vlj;
		REAL ivor_fb_swi = swi_ivor*(6.0 * vlj + shapedelta * (rep * rr23));
		fx = ivor_fb_swi * dx;
		fy = ivor_fb_swi * dy;
		fz = ivor_fb_swi * dz;
	}
}

template<typename REAL>
__host__ __device__
void ChargePotForce(
		REAL const& dr2_inv,
		REAL const& dx,
		REAL const& dy,
		REAL const& dz,
		REAL const& chargeLigRec,
		REAL const& swi,
		Dielec const& dielec,
		REAL& fx,
		REAL& fy,
		REAL& fz,
		REAL& energy)
{
	REAL dd;

	switch(dielec) {
	case Dielec::constant:
		dd = sqrt(dr2_inv) - 1.0 / 50.0;
		break;
	case Dielec::variable:
		dd = dr2_inv - 1.0 / (50.0 * 50.0);
		break;
	}

	/* (cap all distances at 50 A) */
	if (dd < 0)
		dd = 0;

	energy = swi * chargeLigRec * dd;

	switch (dielec) {
	case Dielec::constant:
		if (dd <= 0) {
			fx = 0;
			fy = 0;
			fz = 0;
		} else {
			double et2;
			et2 = swi * chargeLigRec * sqrt(dr2_inv);
			fx = et2 * dx;
			fy = et2 * dy;
			fz = et2 * dz;
		}
		break;
	case Dielec::variable:
		if (dd <= 0) {
			fx = 0;
			fy = 0;
			fz = 0;
		} else {
			double et2;
			et2 = swi * chargeLigRec * dr2_inv;
			fx = 2 * et2 * dx;
			fy = 2 * et2 * dy;
			fz = 2 * et2 * dz;
		}
		break;
	}

}

}  // namespace as



#endif /* SRC_FORCEFIELD_H_ */
