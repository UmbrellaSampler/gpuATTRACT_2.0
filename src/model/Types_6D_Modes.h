/*
 * Types_6D_modes.h
 *
 *  Created on: Nov 4, 2017
 *      Author: glenn
 */

#ifndef SRC_MODEL_TYPES_6D_MODES_H_
#define SRC_MODEL_TYPES_6D_MODES_H_


#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "GenericTypes.h"
#include <vector>
#include <Eigen/Core>
#include "TypeConfiguration.h"

namespace as {
#ifndef __CUDACC__ // ostream is not available in nvcc
template<typename REAL>
struct DOF_6D_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, DOF_6D_Modes<REAL> const&);

template<typename REAL>
struct Result_6D_Modes;

template<typename REAL>
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<REAL> const& args);

#endif



//template<typename REAL>
//struct DOF_6D_Modes {
//	using real_t = typename TypeWrapper<REAL>::real_t;
//	using vec3_t = Vec3<real_t>;
//	vec3_t pos;
//	vec3_t ang;
//	real_t modes[MODES_LIGAND_MAX_NUMBER];
//	unsigned int numModes;
//};

struct Common_Modes {
	id_t gridLigId;
	id_t gridRecId;
	id_t ligId;
	id_t recId;
	id_t tableId;
	id_t paramsId;
};

//template<typename REAL>
//class Result_6D_Modes {
//	using real_t = typename TypeWrapper<REAL>::real_t;
//	using vec3_t = Vec3<real_t>;
//	real_t E;
//	vec3_t pos;
//	vec3_t ang;
//	real_t modes[MODES_LIGAND_MAX_NUMBER];
//	unsigned int numModes;
//};


/**
 * brief: These types resolve a problem with the types converter in meta.h. Since the modes of DOFs might be varying
 * (and in multibodydocking even the number of proteins) these types allow to safe their configuration in a simple struct
 * which is the same for every dof.
 * This way the solverbase can operate independently of the uses dof type and convert them safely back from Eigen::Vector to the specific DOF
 * by just setting the DOF configuration.
 */

template<typename REAL>
class DOF_6D_Modes {
public:
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	vec3_t pos;
	vec3_t ang;
	real_t modes[MODES_LIGAND_MAX_NUMBER];
	unsigned int numModes;


	TypeConfiguration getConfiguration(){
		TypeConfiguration config;
		config._dofSize=1;
		config._numModes=new int[1];
		config._numModes[0]=numModes;
		return config;
	}

	Vector getVectorfromDOF(){
		Vector vec;
		vec  << ang.x, ang.y, ang.z,
				pos.x, pos.y , pos.z;
		for(int mode=0;mode< numModes; mode++){vec  << modes[mode];}
		return vec;
	}

	void setDOFfromVector(Vector vec,TypeConfiguration config){
		ang.x=vec(0);
		ang.y=vec(1);
		ang.z=vec(2);
		pos.x=vec(3);
		pos.y=vec(4);
		pos.z=vec(5);
		for (int mode=0; mode<config._numModes[0];mode++){	modes[mode]=vec(6+mode);}
	}
};



template<typename REAL>
class Result_6D_Modes {
public:
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
	real_t E;
	vec3_t pos;
	vec3_t ang;
	real_t modes[MODES_LIGAND_MAX_NUMBER];
	unsigned int numModes;


	TypeConfiguration getConfiguration(){
			TypeConfiguration config;
			config._dofSize=1;
			config._numModes=new int[1];
			config._numModes[0]=numModes;
			return config;
		}

	ObjGrad getVectorfromDOF(){
		ObjGrad vec;
		vec.obj=E;
		vec.grad  << ang.x, ang.y, ang.z,
				pos.x, pos.y , pos.z;
		for(int mode=0;mode< numModes; mode++){vec.grad  << modes[mode];}
		return vec;}

	void setDOFfromVector(ObjGrad vec,TypeConfiguration config){
		E=vec.obj;
		ang.x=vec.grad(0);
		ang.y=vec.grad(1);
		ang.z=vec.grad(2);
		pos.x=vec.grad(3);
		pos.y=vec.grad(4);
		pos.z=vec.grad(5);
		for (int mode=0; mode<config._numModes[0];mode++){	modes[mode]=vec.grad(6+mode);}
	}
};


template<typename REAL>
class DOF_Vector_6D_Modes {
public:
	std::vector<DOF_6D_Modes<REAL>> dof;

	void setSize(int size){dof.reserve(size);}
	TypeConfiguration getConfiguration(){
		TypeConfiguration config;
		config._dofSize=dof.size();
		config._numModes=new int[dof.size()];
		for(int i=0;i<dof.size();i++){ config._numModes[i]=dof[i].numModes; }
		return config;
	}
	Vector getVectorfromDOF(){
		Vector vec;
		for(int i=0; i<dof.size();i++){
			vec << dof[i].ang.x, dof[i].ang.y, dof[i].ang.z,
							dof[i].pos.x, dof[i].pos.y , dof[i].pos.z;
			for(int mode=0;mode< dof[i].numModes; mode++){vec  << dof[i].modes[mode];}
		}
		return vec;
	}
	void setDOFfromVector(Vector vec,TypeConfiguration config){
		int dofcounter=0;
		dof.reserve(config._dofSize());
		for(int i=0; i<config._dofSize();i++){
			dof[i].pos.x=vec(0+dofcounter);
			dof[i].pos.y=vec(1+dofcounter);
			dof[i].pos.z=vec(2+dofcounter);
			dof[i].ang.x=vec(3+dofcounter);
			dof[i].ang.y=vec(4+dofcounter);
			dof[i].ang.z=vec(5+dofcounter);

			for (int mode=0; mode<config._numModes[i];mode++){	dof[i].modes[mode]=vec(6+mode+dofcounter);
			dofcounter+=config._numModes[i];}
			dofcounter+=6;
		}
	}

};

template<typename REAL>
class Result_Vector_6D_Modes {
public:
	std::vector<Result_6D_Modes<REAL>> dof;

	TypeConfiguration getConfiguration(){
		TypeConfiguration config;
		config._dofSize=dof.size();
		config._numModes=new int[dof.size()];
		for(int i=0;i<dof.size();i++){
			config._numModes[i]=dof[i].numModes;
		}
		return config;
	}

	ObjGrad getVectorfromDOF(){
		ObjGrad vec;
		vec.obj=dof.E;
		for(int i=0; i<dof.size();i++){
			vec.grad  << dof[i].ang.x, dof[i].ang.y, dof[i].ang.z,
							dof[i].pos.x, dof[i].pos.y , dof[i].pos.z;
			for(int mode=0;mode< dof[i].numModes; mode++){vec.grad  << dof[i].modes[mode];}
		}
		return vec;
	}

	void setDOFfromVector(ObjGrad vec,TypeConfiguration config){
		int dofcounter=0;
		dof.reserve(config._dofSize());
		for(int i=0; i<config._dofSize();i++){
			dof[i].E=vec.obj;
			dof[i].pos.x=vec.grad(0+dofcounter);
			dof[i].pos.y=vec.grad(1+dofcounter);
			dof[i].pos.z=vec.grad(2+dofcounter);
			dof[i].ang.x=vec.grad(3+dofcounter);
			dof[i].ang.y=vec.grad(4+dofcounter);
			dof[i].ang.z=vec.grad(5+dofcounter);

			for (int mode=0; mode<config._numModes[i];mode++){	dof[i].modes[mode]=vec.grad(6+mode+dofcounter);
			dofcounter+=config._numModes[i];}
			dofcounter+=6;
		}
	}

};






template<typename REAL>
using Types_6D_Modes = GenericTypes<DOF_Vector_6D_Modes<REAL>, Common_Modes, Result_Vector_6D_Modes<REAL>>;

}  // namespace as




#endif /* TYPES_6D_MODES_H_ */
