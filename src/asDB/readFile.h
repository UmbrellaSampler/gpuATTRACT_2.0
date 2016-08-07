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


///*
// ** @brief: read the number of atoms of a protein from a pdb-file.
// */
//unsigned readProteinSizeFromPDB(std::string filename);
//
///*
// ** @brief: deprecated
// ** Reads binary file according to ProteinDesc.
// ** Memory management needs to handled outside
// ** (e.g. by using a shared_ptr or manual deletion).
// */
//as::Protein* createProteinFromDumpFile(std::string filename);
//
//
//
//
///*
// ** @brief: deprecated
// ** Reads binary file according to GridDesc.
// ** Memory management needs to handled outside
// ** (e.g. by using a shared_ptr or manual deletion).
// */
//as::GridUnion* createGridUnionFromDumpFile (std::string filename);
//
///*

//
//
///*
// ** @brief: reads a (ATTRACT-) .dat containing DOFs. The number of elements read
// ** is returned by the parameter numEl. Normal modes are not yet supported.
// ** Memory needs to deleted by the caller. (e.g. by using a shared_ptr or manual deletion)
// ** ToDo: Read normal modes
// */
//void readDOFFromFile(std::string filename, std::vector<std::vector<as::DOF>>& DOF_molecules);
//
//void readEnsembleDOFFromFile(std::string filename, std::vector<std::vector<as::DOF>>& DOF_molecules);
//
///*
// ** @brief: reads the header of a (ATTRACT-) .dat containing DOFs.
// */
//void readDOFHeader(std::string filename, std::vector<asUtils::Vec3f>& pivots,
//		bool& auto_pivot, bool& centered_receptor, bool& centered_ligands);
//
//std::vector<std::string> readFileNamesFromEnsembleList(std::string filename);
//
//std::vector<unsigned> readGridAlphabetFromFile(std::string filename);

} // namespace as









#endif /* ASDB_READ_FILE_H_ */
