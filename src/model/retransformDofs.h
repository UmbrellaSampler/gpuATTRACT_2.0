/*
 * retransformDofs.h
 *
 *  Created on: May 7, 2018
 *      Author: glenn
 */

#ifndef RETRANSFORMDOFS_H_
#define RETRANSFORMDOFS_H_

#include "retransformDofs.tpp"

namespace as{


template
void transformDOF<double>(	bool input_centered_receptor, bool input_centered_ligand,
					bool output_centered_receptor, bool output_centered_ligand,
					Vec3<double> pivot_receptor,Vec3<double> pivot_ligand,
					std::vector<DOF_6D<double>> &dofs);
template
void transformDOF<float>(	bool input_centered_receptor, bool input_centered_ligand,
					bool output_centered_receptor, bool output_centered_ligand,
					Vec3<double> pivot_receptor,Vec3<double> pivot_ligand,
					std::vector<DOF_6D<float>> &dofs);
template
void transformDOF<double>(	bool input_centered_receptor, bool input_centered_ligand,
					bool output_centered_receptor, bool output_centered_ligand,
					Vec3<double> pivot_receptor,Vec3<double> pivot_ligand,
					std::vector<DOF_6D_Modes<double>> &dofs);
template
void transformDOF<float>(	bool input_centered_receptor, bool input_centered_ligand,
					bool output_centered_receptor, bool output_centered_ligand,
					Vec3<double> pivot_receptor,Vec3<double> pivot_ligand,
					std::vector<DOF_6D_Modes<float>> &dofs);
}

#endif /* RETRANSFORMDOFS_H_ */
