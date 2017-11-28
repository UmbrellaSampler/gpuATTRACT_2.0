/*
 * TypeConfiguration.h
 *
 *  Created on: Nov 24, 2017
 *      Author: glenn
 */

#ifndef TYPECONFIGURATION_H_
#define TYPECONFIGURATION_H_

#include <Eigen/Core>



namespace as {

const unsigned MODES_LIGAND_MAX_NUMBER=10;
const unsigned MODES_RECEPTOR_MAX_NUMBER=10;

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Scalar = Eigen::VectorXd::Scalar;

struct TypeConfiguration{
	int _dofSize;
	int* _numModes;
};


struct ObjGrad {
	double obj; // function value
	Vector grad; // gradients
};
/**
 * @brief: this class was mainly build for the Requesthandler which also converts from Eigen::Vector to the dofs.
 * informations like number of Proteins or number of modes are not contained in this type.
 * In this case an object of typeconfiguration can be called, which saves all relevant information and upon calling returns Input or outputtypes
 */

//
//
//class TypeConfiguration{
//private:
//	int _dofSize;
//	int* _numModes;
//	bool hasModes;
//	bool isVector;
//public:
////	template<typename INPUT >
////	void setConfiguration(INPUT input);
//
//	template<typename RETURNTYPE>
//	RETURNTYPE createResultDOFfromConfiguration();
//
//	template<typename RETURNTYPE>
//		RETURNTYPE createDOFfromVector(Vector vec);
//
//
//
//void setConfiguration(std::vector<DOF_6D_Modes<float>> input){
//	_numModes=new int[input.size()];
//	_dofSize=input.size();
//	hasModes=true;
//	isVector=true;
//			for(int i=0;i<input.size();i++){
//				_numModes[i]=input.at(i).numModes;
//			}
//	}
//
//	void setConfiguration(std::vector<DOF_6D<float>> input){
//		_numModes=new int[input.size()];
//		_dofSize=input.size();
//		hasModes=true;
//		isVector=true;
//
//		}
//
//	void setConfiguration(const DOF_6D<float> &input){
//		_dofSize=1;
//		_numModes=NULL;
//		hasModes=false;
//		isVector=false;
//	}
//
//
//	void setConfiguration(const DOF_6D<double> &input){
//		_dofSize=1;
//		_numModes=NULL;
//		hasModes=false;
//		isVector=false;
//	}
//
//	void setConfiguration(const Result_6D<float> &input){
//		_dofSize=1;
//		_numModes=NULL;
//		hasModes=false;
//		isVector=false;
//	}
//
//
//	void setConfiguration(const Result_6D<double> &input){
//		_dofSize=1;
//		_numModes=NULL;
//		hasModes=false;
//		isVector=false;
//	}
//
//	void setConfiguration(DOF_6D_Modes<float> input){
//			_dofSize=1;
//			_numModes=new int[1];
//			_numModes=NULL;
//			hasModes=true;
//			_numModes[0]=input.numModes;
//			isVector=false;
//		}
//};
//	template<typename RETURNTYPE>
//	RETURNTYPE TypeConfiguration::createResultDOFfromConfiguration(){
//		if(isVector){RETURNTYPE dof(_dofSize);
//		if(hasModes){
//			if(isVector){
//			for(int i=0;i<_dofSize;i++){
//				dof[i].numModes=_numModes[i];
//			}}
//			else{dof.numModes=_numModes[0];}
//
//			return dof;
//		}
//		else{return dof;}
//		}
//		else{
//		RETURNTYPE dof;
//
//		if(hasModes){
//			dof.numModes=_numModes[0];}
//		else{dof.numModes=_numModes[0];}
//		return dof;}
//	}
//
//	/////////////////////	CREATE VECTOR
//	template<typename RETURNTYPE>
//	RETURNTYPE TypeConfiguration::createDOFfromVector(Vector vec){
////		if(isVector){
////			RETURNTYPE dof(_dofSize);
////
////			int dofcounter=0;
////			for(int i=0;i<_dofSize;i++){
////
////				dof[i].ang.x=vec(0+dofcounter);
////				dof[i].ang.y=vec(1+dofcounter);
////				dof[i].ang.z=vec(2+dofcounter);
////				dof[i].pos.x=vec(3+dofcounter);
////				dof[i].pos.y=vec(4+dofcounter);
////				dof[i].pos.z=vec(5+dofcounter);
////				if(hasModes){
////					for (int mode=0; mode<_numModes[i];mode++){	dof[i].modes[mode]=vec(6+mode+dofcounter);}
////					dofcounter+_numModes[i];
////				}
////			dofcounter+=6;
////		}
////			return dof;
////		}
////		else{
//			RETURNTYPE dof;
//
//
//		dof.ang.x=vec(0);
//		dof.ang.y=vec(1);
//		dof.ang.z=vec(2);
//		dof.pos.x=vec(3);
//		dof.pos.y=vec(4);
//		dof.pos.z=vec(5);
//		if(hasModes){
//			for (int mode=0; mode<_numModes[0];mode++){	dof.modes[mode]=vec(6+mode);}
//		}
//		return dof;
//	//	}
//
//			}
//
//
//















/////////NO  VECTOR ///////  MODES
//
//	template<>
//	void TypeConfiguration::setConfiguration(DOF_6D_Modes<float> input){
//	_numModes=new int[1];
//	_numModes[0]=input.numModes;
//	}
//	/////////////////////	CREATE DOF
//	template<>
//	DOF_6D_Modes<float> TypeConfiguration::createInputDOFfromConfiguration(){
//			DOF_6D_Modes<float> dof;
//			dof.numModes=_numModes[0];
//			return dof;
//		}
//
//	template<>
//	Result_6D_Modes<float> TypeConfiguration::createResultDOFfromConfiguration(){
//		Result_6D_Modes<float> dof;
//		dof.numModes=_numModes[0];
//		return dof;
//		}
//
//
//	/////////////////////	CREATE VECTOR
//	template<>
//	Result_6D_Modes<float> TypeConfiguration::createResultDOFfromVector(Vector vec){
//		Result_6D_Modes<float> dof;
//			dof.ang.x=vec(0);
//			dof.ang.y=vec(1);
//			dof.ang.z=vec(2);
//			dof.pos.x=vec(3);
//			dof.pos.y=vec(4);
//			dof.pos.z=vec(5);
//			for (int mode=0; mode<_numModes[0];mode++){	dof.modes[mode]=vec(6+mode);}
//		return dof;
//	}
/////////NO  VECTOR /////// NO MODES
//
//
//	template<>
//		void TypeConfiguration::setConfiguration(DOF_6D<float> input){
//		_dofSize=1;
//		_numModes=new int[1];
//		_numModes[0]=0;
//		}
//
//
//	template<>
//	DOF_6D<float> TypeConfiguration::createInputDOFfromConfiguration(){
//			DOF_6D<float> dof;
//			return dof;
//		}
//
//	template<>
//	Result_6D<float> TypeConfiguration::createResultDOFfromConfiguration(){
//			Result_6D<float> dof;
//			return dof;
//		}
//
//
//	/////////////////////	CREATE VECTOR
//	template<>
//	Result_6D<float> TypeConfiguration::createResultDOFfromVector(Vector vec){
//		Result_6D<float> dof;
//			dof.ang.x=vec(0);
//			dof.ang.y=vec(1);
//			dof.ang.z=vec(2);
//			dof.pos.x=vec(3);
//			dof.pos.y=vec(4);
//			dof.pos.z=vec(5);
//		return dof;
//	}










//		}
//template<>
//TypeConfiguration::setConfiguration(std::vector<Result_6D_Modes<double>> input):_dofSize(input.size()),_numModes(new int[input.size()]){
//			for(int i=0;i<input.size();i++){
//				_numModes[i]=input.at(i).numModes;
//			}
//		}
//
//
//
//
//
//
//
//template<>
//		std::vector<DOF_6D<double>> TypeConfiguration::createInputDOFfromConfiguration(){
//					std::vector<DOF_6D<double>> dof(_dofSize);
//					return dof;
//				}
//template<>
//		std::vector<Result_6D<double>> TypeConfiguration::createResultDOFfromConfiguration(){
//					std::vector<Result_6D<double>> dof(_dofSize);
//					return dof;
//				}
//
//
//
//
//
//
//
//
//
//
//
//
//template<>
//		TypeConfiguration::TypeConfiguration(DOF_6D<double> input){}
//template<>
//		TypeConfiguration::TypeConfiguration(Result_6D<double> input){}
//template<>
//			DOF_6D<double> TypeConfiguration::createInputDOFfromConfiguration(){
//				DOF_6D<double> dof;
//				return dof;
//			}
//template<>
//			Result_6D<double> TypeConfiguration::createResultDOFfromConfiguration(){
//				Result_6D<double> dof;
//				return dof;
//				}
//
//
//
//
//
//template<>
//			TypeConfiguration::TypeConfiguration(DOF_6D_Modes<double> input):_numModes(new int[1]){
//	_numModes
//}
//
//template<>
//			DOF_6D_Modes<double> TypeConfiguration::createInputDOFfromConfiguration(){
//						DOF_6D_Modes<double> dof;
//						dof.numModes=_numModes;
//						return dof;
//					}
//
//template<>
//			TypeConfiguration::TypeConfiguration(Result_6D_Modes<double> input):_numModes(input.numModes){}
//
//template<>
//			DOF_6D_Modes<double> TypeConfiguration::createInputDOFfromConfiguration(){
//				DOF_6D_Modes<double> dof;
//				dof.numModes=_numModes;
//				return dof;
//			}
//template<>
//			Result_6D_Modes<double> TypeConfiguration::createResultDOFfromConfiguration(){
//				Result_6D_Modes<double> dof;
//				dof.numModes=_numModes;
//				return dof;
//				}
//
//
//
//
//
//
//template<typename INPUT, typename RETURNTYPE>
//class TypeConfiguration{
//private:
//	int _dofSize;
//	int* _numModes;
//public:
//	TypeConfiguration(INPUT input);
//	RETURNTYPE createDOFfromConfiguration();
//};
//
//
//
////template <>
////class TypeConfiguration<std::vector<DOF_6D_Modes<float>> , std::vector<Result_6D_Modes<float>> >;
////template <>
////class TypeConfiguration<std::vector<DOF_6D_Modes<double>> , std::vector<Result_6D_Modes<double>> >;
////template <>
////class TypeConfiguration<std::vector<DOF_6D<float>> , std::vector<Result_6D<float>> >;
////template <>
////class TypeConfiguration<std::vector<DOF_6D<double>> , std::vector<Result_6D<double>> >;
////
////template <>
////class TypeConfiguration<DOF_6D_Modes<float> , Result_6D_Modes<float> >;
////template <>
////class TypeConfiguration<DOF_6D_Modes<double> , Result_6D_Modes<double> >;
////template <>
////class TypeConfiguration<DOF_6D<float>, Result_6D<float>>;
////template <>
////class TypeConfiguration<DOF_6D<double>, Result_6D<double>>;
//
//
//template <>
//class TypeConfiguration<std::vector<DOF_6D_Modes<double>> , std::vector<Result_6D_Modes<double>> >{
//	private:
//		int _dofSize;
//		int* _numModes;
//
//	public:
//	TypeConfiguration(std::vector<DOF_6D_Modes<double>> input):_dofSize(input.size()),_numModes(new int[input.size()]){
//			for(int i=0;i<input.size();i++){
//				_numModes[i]=input.at(i).numModes;
//			}
//		}
//	TypeConfiguration(std::vector<Result_6D_Modes<double>> input):_dofSize(input.size()),_numModes(new int[input.size()]){
//			for(int i=0;i<input.size();i++){
//				_numModes[i]=input.at(i).numModes;
//			}
//		}
//
//	std::vector<DOF_6D_Modes<double>> createInputDOFfromConfiguration(){
//			std::vector<DOF_6D_Modes<double>> dof(_dofSize);
//			for(int i=0;i<_dofSize;i++){
//				dof[i].numModes=_numModes[i];
//			}
//			return dof;
//		}
//
//	std::vector<Result_6D_Modes<double>> createResultDOFfromConfiguration(){
//			std::vector<Result_6D_Modes<double>> dof(_dofSize);
//			for(int i=0;i<_dofSize;i++){
//				dof[i].numModes=_numModes[i];
//			}
//			return dof;
//		}
//	};
//
//
//
//template <>
//class TypeConfiguration<std::vector<DOF_6D<double> >, std::vector<Result_6D<double>> >{
//private:
//	int _dofSize;
//	int* _numModes;
//public:
//	TypeConfiguration(std::vector<DOF_6D<float>> input){
//		_dofSize=input.size();
//		_numModes=0;
//		}
//	TypeConfiguration(std::vector<DOF_6D<double>> input){
//		_dofSize=input.size();
//		_numModes=0;
//		}
//	std::vector<DOF_6D<double>> createInputDOFfromConfiguration(){
//				std::vector<DOF_6D<double>> dof(_dofSize);
//				return dof;
//			}
//
//	std::vector<Result_6D<double>> createResultDOFfromConfiguration(){
//				std::vector<Result_6D<double>> dof(_dofSize);
//				return dof;
//			}
//};
//
//
//template <>
//class TypeConfiguration<DOF_6D<double> , Result_6D<double> >{
//public:
//	TypeConfiguration(DOF_6D<double> input){}
//	TypeConfiguration(Result_6D<double> input){}
//
//	DOF_6D<double> createInputDOFfromConfiguration(){
//		DOF_6D<double> dof;
//		return dof;
//	}
//
//	Result_6D<double> createResultDOFfromConfiguration(){
//		Result_6D<double> dof;
//		return dof;
//		}
//};
//
//
//template <>
//class TypeConfiguration<DOF_6D_Modes<double> , Result_6D_Modes<double> >{
//private:
//	int _numModes;
//public:
//	TypeConfiguration(DOF_6D_Modes<double> input):_numModes(input.numModes){}
//	TypeConfiguration(Result_6D_Modes<double> input):_numModes(input.numModes){}
//
//	DOF_6D_Modes<double> createInputDOFfromConfiguration(){
//		DOF_6D_Modes<double> dof;
//		dof.numModes=_numModes;
//		return dof;
//	}
//
//	Result_6D_Modes<double> createResultDOFfromConfiguration(){
//		Result_6D_Modes<double> dof;
//		dof.numModes=_numModes;
//		return dof;
//		}
//};
//
//
//template <>
//class TypeConfiguration<std::vector<DOF_6D_Modes<float>> , std::vector<Result_6D_Modes<float>> >{
//	private:
//		int _dofSize;
//		int* _numModes;
//
//	public:
//	TypeConfiguration(std::vector<DOF_6D_Modes<float>> input):_dofSize(input.size()),_numModes(new int[input.size()]){
//			for(int i=0;i<input.size();i++){
//				_numModes[i]=input.at(i).numModes;
//			}
//		}
//	TypeConfiguration(std::vector<Result_6D_Modes<float>> input):_dofSize(input.size()),_numModes(new int[input.size()]){
//			for(int i=0;i<input.size();i++){
//				_numModes[i]=input.at(i).numModes;
//			}
//		}
//
//	std::vector<DOF_6D_Modes<float>> createInputDOFfromConfiguration(){
//			std::vector<DOF_6D_Modes<float>> dof(_dofSize);
//			for(int i=0;i<_dofSize;i++){
//				dof[i].numModes=_numModes[i];
//			}
//			return dof;
//		}
//
//	std::vector<Result_6D_Modes<float>> createResultDOFfromConfiguration(){
//			std::vector<Result_6D_Modes<float>> dof(_dofSize);
//			for(int i=0;i<_dofSize;i++){
//				dof[i].numModes=_numModes[i];
//			}
//			return dof;
//		}
//	};
//
//
//
//template <>
//class TypeConfiguration<std::vector<DOF_6D<float> >, std::vector<Result_6D<float>> >{
//private:
//	int _dofSize;
//	int* _numModes;
//public:
//	TypeConfiguration(std::vector<DOF_6D<float>> input){
//		_dofSize=input.size();
//		_numModes=0;
//		}
//	TypeConfiguration(std::vector<DOF_6D<double>> input){
//		_dofSize=input.size();
//		_numModes=0;
//		}
//	std::vector<DOF_6D<float>> createInputDOFfromConfiguration(){
//				std::vector<DOF_6D<float>> dof(_dofSize);
//				return dof;
//			}
//
//	std::vector<Result_6D<float>> createResultDOFfromConfiguration(){
//				std::vector<Result_6D<float>> dof(_dofSize);
//				return dof;
//			}
//};
//
//
//template <>
//class TypeConfiguration<DOF_6D<float> , Result_6D<float> >{
//public:
//	TypeConfiguration(DOF_6D<float> input){}
//	TypeConfiguration(Result_6D<float> input){}
//
//	DOF_6D<float> createInputDOFfromConfiguration(){
//		DOF_6D<float> dof;
//		return dof;
//	}
//
//	Result_6D<float> createResultDOFfromConfiguration(){
//		Result_6D<float> dof;
//		return dof;
//		}
//};
//
//
//template <>
//class TypeConfiguration<DOF_6D_Modes<float> , Result_6D_Modes<float> >{
//private:
//	int _numModes;
//public:
//	TypeConfiguration(DOF_6D_Modes<float> input):_numModes(input.numModes){}
//	TypeConfiguration(Result_6D_Modes<float> input):_numModes(input.numModes){}
//
//	DOF_6D_Modes<float> createInputDOFfromConfiguration(){
//		DOF_6D_Modes<float> dof;
//		dof.numModes=_numModes;
//		return dof;
//	}
//
//	Result_6D_Modes<float> createResultDOFfromConfiguration(){
//		Result_6D_Modes<float> dof;
//		dof.numModes=_numModes;
//		return dof;
//		}
//};






}// end namespace as

#endif /* TYPECONFIGURATION_H_ */
