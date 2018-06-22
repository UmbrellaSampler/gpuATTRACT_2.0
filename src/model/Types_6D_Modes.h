/*
 * Types_6D_modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TYPES_6D_MODES_H_
#define SRC_MODEL_TYPES_6D_MODES_H_

#include "Types_6D_Config.h"
#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "GenericTypes.h"
#include "Types_6D.h"

namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_6D_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_6D_Modes<REAL> const& args);

template<typename REAL>
struct Result_6D_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<REAL> const& args);

#endif

template<typename REAL>
struct DOF_6D_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	DOF_6D<real_t> _6D;
	real_t modesRec[MODES_MAX_NUMBER];
	real_t modesLig[MODES_MAX_NUMBER];
	vec3_t get_pos(){
			return (_6D.pos);
		}
		void set_pos(real_t x,real_t y,real_t z ){
			_6D.pos.x = x;
			_6D.pos.y = y;
			_6D.pos.z = z;
					}

};

struct Common_Modes {
	id_t gridIdRec;
	id_t gridIdLig;
	id_t ligId;
	id_t recId;
	id_t tableId;
	id_t paramsId;
	static unsigned int numModesRec;
	static unsigned int numModesLig;
	Vec3<double> pivotRec;
	Vec3<double> pivotLig;
	bool centeredLig;
	bool centeredRec;
	double radius_cutoff;
};

template<typename REAL>
struct Result_6D_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	Result_6D<real_t> _6D;
	real_t modesRec[MODES_MAX_NUMBER];
	real_t modesLig[MODES_MAX_NUMBER];
	// This method is just for testing purposes in emattract
	REAL get_Energy(){ return _6D.E;}
};

template<typename REAL>
using Types_6D_Modes = GenericTypes<DOF_6D_Modes<REAL>, Common_Modes, Result_6D_Modes<REAL>>;

}  // namespace as




#endif /* TYPES_6D_MODES_H_ */
