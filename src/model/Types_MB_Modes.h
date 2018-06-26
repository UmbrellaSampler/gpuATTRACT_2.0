/*
 * Types_MB_modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TYPES_MB_MODES_H_
#define SRC_MODEL_TYPES_MB_MODES_H_

#include "Types_6D_Config.h"
#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "GenericTypes.h"
#include <map>
#include <vector>


namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_MB_Modes;

template<typename REAL>
struct DOF_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_Modes<REAL> const& args);

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_MB_Modes<REAL> const& args);

template<typename REAL>
struct Result_MB_Modes;



template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_MB_Modes<REAL> const& args);

#endif


template<typename REAL>
struct DOF_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	vec3_t pos;
	vec3_t ang;
	real_t modes[MODES_MAX_NUMBER];
};

template<typename REAL>
struct DOF_MB_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	DOF_Modes<real_t> protein[NUM_MAX_PROTEIN];
	vec3_t get_pos(unsigned idx_protein ){
		return (protein[idx_protein].pos);
	}
	void set_pos(real_t x,real_t y,real_t z, unsigned idx_protein  ){
		protein[idx_protein].pos.x = x;
		protein[idx_protein].pos.y = y;
		protein[idx_protein].pos.z = z;
	}

};

struct ProtConfig{
	using real_t = double;
	using vec3_t = Vec3<double>;
	id_t gridId;
	id_t proteinId;
	unsigned idxModes;
	bool centered;
	vec3_t pivot;
	template< typename REAL>
	ProtConfig( id_t id_grid, id_t id_protein, id_t idx_mode, bool center, Vec3<REAL> piv): gridId(id_grid), proteinId(id_protein), idxModes(idx_mode), centered(center), pivot(piv.x,piv.y,piv.z){}
	//template <typename REAL>
//	ProtConfig( id_t id_grid, id_t id_protein, id_t idx_mode, bool center, Vec3<REAL> piv){
//		gridId = id_grid;
//		proteinId = id_protein;
//		idxModes = idx_mode;
//		centered = center;
//		pivot = vec3_t( piv.x,piv.y,piv.z);
//	}
};
struct Common_MB_Modes {
	id_t tableId;
	id_t paramsId;
	static unsigned int numModes[NUM_MAX_PROTEIN];
	double radius_cutoff;
	//std::map<unsigned,ProtConfig> proteins;
	std::vector<ProtConfig> proteins;
	static unsigned int numProteins;
	Vec3<double> getPivot( unsigned idx_protein)
	{
		return proteins[idx_protein].pivot;
	}

};

template<typename REAL>
struct Result_MB_Modes {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	DOF_Modes<real_t> protein[NUM_MAX_PROTEIN];
	real_t E;
	REAL get_Energy(){ return E;}

};

template<typename REAL>
using Types_MB_Modes = GenericTypes<DOF_MB_Modes<REAL>, Common_MB_Modes, Result_MB_Modes<REAL>>;

}  // namespace as




#endif /* TYPES_MB_MODES_H_ */
