/*
 * debug_functions.h
 *
 *  Created on: Jun 25, 2018
 *      Author: glenn
 */

#ifndef DEBUG_FUNCTIONS_H_
#define DEBUG_FUNCTIONS_H_

#include <sstream>
#include "Types_6D_Modes.h"
#include "Types_6D.h"
#include <iomanip>

namespace as
{


template <typename REAL,typename DOF_T>
std::string getDebugPath(bool gpu )
{
	std::string path_debug = "/home/glenn/Documents/test_attract/output";
	std::string path_out = path_debug;
	if(gpu){path_out += "/GPU";}
	else{path_out += "/CPU";}

	if (  (std::is_same<DOF_T, Result_6D_Modes<REAL>>::value || std::is_same<DOF_T, DOF_6D_Modes<REAL>>::value)){
		path_out += "_Modes";}

	return path_out;

}

template <typename REAL,typename DOF_T>
std::string getDebugFile(bool gpu )
{
	std::string file_out = getDebugPath< REAL, DOF_T>(gpu );
	if (  (std::is_same<DOF_T, Result_6D_Modes<REAL>>::value || std::is_same<DOF_T, Result_6D<REAL>>::value)){
		file_out += "/result.dat";
	}
	else {file_out += "/dof.dat";	}
	return file_out;

}


template< typename REAL , typename DOF_T>
void print_header( std::stringstream& os,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	os <<" " << "rec_rot_x"<< " " << "rec_rot_y"<< " " << "rec_rot_z"<< " " << "rec_trans_x"<< " " << "rec_trans_y"<< " " << "rec_trans_z" ;
	for ( int i = 0; i < Common_Modes::numModesRec; ++i)
	{
		os << " " << "rec_mode_"<< i+1;
	}
	os << " " << "lig_rot_x"<< " " << "lig_rot_y"<< " " << "lig_rot_z"<< " " << "lig_trans_x"<< " " << "lig_trans_y"<< " " << "lig_trans_z";
	for ( int i = 0; i < Common_Modes::numModesLig; ++i)
	{
		os << " " << "lig_mode_"<< i+1;
	}
}

template< typename REAL , typename DOF_T>
void print_header( std::stringstream& os,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0 )
{
	os << " " << "lig_rot_x"<< " " << "lig_rot_y"<< " " << "lig_rot_z"<< " " << "lig_trans_x"<< " " << "lig_trans_y"<< " " << "lig_trans_z";
}


template< typename REAL , typename DOF_T>
void print_header( std::stringstream& os,
		typename std::enable_if<std::is_same< DOF_T, Result_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	os <<" "<< "energy_total";
	print_header<REAL, DOF_6D_Modes<REAL>>( os);

}

template< typename REAL , typename DOF_T>
void print_header( std::stringstream& os,
		typename std::enable_if<std::is_same< DOF_T, Result_6D<REAL> >::value, void>::type* dummy = 0 )
{
	os <<" "<< "energy_total";
	print_header<REAL, DOF_6D<REAL>>( os);
}

template< typename REAL , typename DOF_T>
void print_energy( std::stringstream& os, DOF_T res,
		typename std::enable_if<std::is_same< DOF_T, Result_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	auto w = std::setw(20);
	os  << w << res._6D.E;
}

template< typename REAL , typename DOF_T>
void print_energy( std::stringstream& os, DOF_T const res,
		typename std::enable_if<std::is_same< DOF_T, Result_6D<REAL> >::value, void>::type* dummy = 0 )
{
	auto w = std::setw(20);
	os  << w << res.E;

}

template< typename REAL , typename DOF_T>
void print_enGrad( std::stringstream& os, DOF_T const  res,
		typename std::enable_if<std::is_same< DOF_T, Result_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	auto w = std::setw(20);
	os  << w << res._6D.E;
	os  << w <<  " " << w;
	os 	<< w << 0 << w << 0 << w << 0;
	os  << w << 0 << w << 0 << w << 0;
	for(int mode=0;mode<Common_Modes::numModesRec;mode++){
		os<< w << res.modesRec[mode];
	}
	os  << w <<  " " << w;

	os	<< w << res._6D.ang.x << w << res._6D.ang.y << w << res._6D.ang.z;
	os	<< w << res._6D.pos.x << w << res._6D.pos.y << w << res._6D.pos.z;
	for(int mode=0;mode<Common_Modes::numModesLig;mode++){
		os<< w << res.modesLig[mode];
	}
	os  << w <<  " " << w;

}

template< typename REAL , typename DOF_T>
void print_enGrad( std::stringstream& os, DOF_T const res,
		typename std::enable_if<std::is_same< DOF_T, Result_6D<REAL> >::value, void>::type* dummy = 0 )
{
	auto w = std::setw(20);
	os  << w << res.E;
	os  << w <<  " " << w;
	os	<< w << res.ang.x << w << res.ang.y << w << res.ang.z;
	os	<< w << res.pos.x << w << res.pos.y << w << res.pos.z;
	os  << w <<  " " << w;

}

template< typename REAL , typename DOF_T>
void print_dof( std::stringstream& os, DOF_T const dof,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D_Modes<REAL> >::value, void>::type* dummy = 0 )
{
	auto w = std::setw(20);
	os 	<< w << 0 << w << 0 << w << 0;
	os  << w << 0 << w << 0 << w << 0;
	for(int mode=0;mode<Common_Modes::numModesRec;mode++){
		os<< w << dof.modesRec[mode];
	}
	os  << w <<  " " << w;

	os	<< w << dof._6D.ang.x << w << dof._6D.ang.y << w << dof._6D.ang.z;
	os	<< w << dof._6D.pos.x << w << dof._6D.pos.y << w << dof._6D.pos.z;
	for(int mode=0;mode<Common_Modes::numModesLig;mode++){
		os<< w << dof.modesLig[mode];
	}
	os  << w <<  " " << w;
}


template< typename REAL , typename DOF_T>
void print_dof( std::stringstream& os, DOF_T const dof,
		typename std::enable_if<std::is_same< DOF_T, DOF_6D<REAL> >::value, void>::type* dummy = 0 )
{
	auto w = std::setw(20);
	os	<< w << dof.ang.x << w << dof.ang.y << w << dof.ang.z;
	os	<< w << dof.pos.x << w << dof.pos.y << w << dof.pos.z;
	os  << w <<  " " << w;
}


}

#endif /* DEBUG_FUNCTIONS_H_ */
