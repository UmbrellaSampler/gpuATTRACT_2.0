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

#ifndef ASDB_READ_FILE_H_
#define ASDB_READ_FILE_H_

#include <string>
#include <memory>
#include <vector>
#include "IOException.h"
#include "nativeTypesWrapper.h"
#include "Vec3.h"
#include "Types_6D.h"

namespace as {

template<typename REAL>
class Protein;

template<typename REAL>
class GridUnion;

template<typename REAL>
class ParamTable;

/*
 ** @brief: Creates Protein object, reads pdb and assignes values.
 */
template<typename REAL>
std::shared_ptr<Protein<REAL>> createProteinFromPDB (std::string filename);

/*
 ** @brief: reads Modes from Modes files
 */
template<typename REAL>
void readHMMode(std::shared_ptr<Protein<REAL>> prot, std::string modeFileName);

/*
 ** @brief: as above but user provides Protein pointer.
 */
template<typename REAL>
void readProteinFromPDB(std::shared_ptr<Protein<REAL>>, std::string filename);

/*
 ** @brief: Creates GridUnion from ATTRACT grid file (original
 ** format).
 */
template<typename REAL>
std::shared_ptr<GridUnion<REAL>> createGridFromGridFile(std::string filename);

/*
 ** @brief: as above but user provides GridUnion pointer.
 */
template<typename REAL>
void readGridFromGridFile(std::shared_ptr<GridUnion<REAL>>, std::string filename);

/*
 ** @brief: reads an ATTRACT forcefield parameter file and creates an
 ** ATTRACT parameter table object.
 */
template<typename REAL>
std::shared_ptr<ParamTable<REAL>> createParamTableFromFile(std::string filename);

/*
 ** @brief: as above but user provides AttrParamTable pointer.
 */
template<typename REAL>
void readParamTableFromFile(std::shared_ptr<ParamTable<REAL>>, std::string filename);

struct AttractEnGrad {
	double E;
	double E_VdW;
	double E_El;
	double3 pos;
	double3 ang;
};

std::vector<AttractEnGrad> readEnGradFromFile(std::string filename);

template<typename REAL>
class DOFHeader {
	using real_t = typename TypeWrapper<REAL>::real_t;
	using vec3_t = Vec3<real_t>;
public:
	std::vector<vec3_t> pivots;
	bool auto_pivot;
	bool centered_receptor;
	bool centered_ligands;
};

constexpr unsigned DOF_MAX_NUMBER = 16;

class DOF {
public:
	DOF_6D<double> _6D;
	double dofs[DOF_MAX_NUMBER-6];
	unsigned numDofs;
};

std::vector<std::vector<DOF>> readDOF(std::string filename);

template<typename REAL>
DOFHeader<REAL> readDOFHeader(std::string filename);

std::vector<std::string> readFileNamesFromEnsembleList(std::string filename);

std::vector<unsigned> readGridAlphabetFromFile(std::string filename);

} // namespace as









#endif /* ASDB_READ_FILE_H_ */
