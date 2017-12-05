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



//#include "grid_orig.h"

#include "../fileIO/readFile.h"

#include <sstream>
#include <fstream>
#include <memory>
#include <cstring>

#include <iostream>


#include <iterator>

#include <vector>
#include <algorithm>

#include "../fileIO/grid_orig.h"
#include "Vec3.h"
#include "Protein.h"
#include "GridUnion.h"
#include "ParamTable.h"


namespace as {

static void errorDOFFormat(std::string filename) {
	std::cerr << "Error reading file " << filename
			<< ". DOF-Format is not supported." << std::endl;
}

/**
 * Converts a line in pdb-file to list of column entries
 *
 * @param line:
 * @return list of pdb column entries
 */
static std::vector<std::string> pdbLine2Strings(std::string const& line) {
	using namespace std;
	istringstream iss(line);
	/* char offset in file: X, Y, Z, type, charge */
	const int idx[][2] = {{30,38},{38,46},{46,54},{54,59},{59,67}};
	const int size = 5;
	vector<string> tokens;

	for (int i = 0; i < size; ++i) {
		int len = idx[i][1]- idx[i][0];
		string str(line.substr(idx[i][0], len));
		// remove white spaces
		str.erase(remove_if(str.begin(), str.end(), ::isspace),str.end());
		tokens.push_back(str);
	}
	return tokens;
}

static std::vector<std::string> line2Strings(std::string const& line) {
	using namespace std;
	istringstream iss(line);
	return vector<string> { istream_iterator<string> { iss },
					istream_iterator<string> { } };
}

template<typename REAL>
std::shared_ptr<Protein<REAL>> createProteinFromPDB(std::string filename) {
	std::shared_ptr<Protein<REAL>> prot = std::make_shared<Protein<REAL>>();
	readProteinFromPDB(prot, filename);
	return prot;

}

template<typename REAL>
void readProteinFromPDB(std::shared_ptr<Protein<REAL>> prot, std::string filename) {
	using namespace std;

	ifstream file(filename);
	string line;

	vector<REAL> posX;
	vector<REAL> posY;
	vector<REAL> posZ;
	vector<unsigned> types;
	vector<REAL> charges;

	if (file.is_open()) {
		getline(file, line);
		while (!file.eof()) {

			if (line.empty()) {
				getline(file, line);
				continue;
			}

			/* split line to vector of strings */
			vector<string> tokens = pdbLine2Strings(line);

//			for(auto str : tokens) {
//				cout << str << " ";
//			}
//			cout << endl;

			REAL x, y, z;
			unsigned type;
			REAL charge;

			stringstream(tokens[0]) >> x;
			stringstream(tokens[1]) >> y;
			stringstream(tokens[2]) >> z;
			stringstream(tokens[3]) >> type;
			stringstream(tokens[4]) >> charge;

			/* store values in container */
			posX.push_back(x);
			posY.push_back(y);
			posZ.push_back(z);
			types.push_back(type);
			charges.push_back(charge);

			getline(file, line);

		}
	} else {
		cerr << "Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	/* concatenate pos-vectors */
	posX.insert(posX.end(), posY.begin(), posY.end());
	posX.insert(posX.end(), posZ.begin(), posZ.end());

	/* copy values to protein buffers */
	prot->setNumAtoms(types.size()); // set the number of atoms
	REAL* protbufPos = prot->getOrCreatePosPtr();
	std::copy(posX.begin(), posX.end(), protbufPos);
	unsigned* protbufType = prot->getOrCreateTypePtr();
	std::copy(types.begin(), types.end(), protbufType);
	REAL* protbufCharge = prot->getOrCreateChargePtr();
	std::copy(charges.begin(), charges.end(), protbufCharge);

	/* set tag */
	prot->setTag(filename);

	file.close();
}

template<typename REAL>
std::shared_ptr<GridUnion<REAL>> createGridFromGridFile(std::string filename) {
	std::shared_ptr<GridUnion<REAL>> grid = std::make_shared<GridUnion<REAL>>();
	readGridFromGridFile(grid, filename);
	return grid;
}

template<typename REAL>
void readHMMode(std::shared_ptr<Protein<REAL>> prot, std::string modeFileName) {
	using namespace std;
	int numModes= prot->numModes();
	int numAtoms=prot->numAtoms();

	vector<REAL> modeX;	vector<REAL> modeY;	vector<REAL> modeZ;
	REAL modeVal[numModes];
	REAL eigVal[numModes];
	REAL *protModes=NULL;
	protModes = new REAL[3*numModes*numAtoms];

	bool changeMode=true;
	bool isData=false;
	int idxPos=0;
	int idxLine=0;
	int modeIdx=0;

	vector<std::string> tokens;
	string line;
	ifstream modeFile(modeFileName);

	if (!modeFile.is_open()){	perror(("error while opening file " + modeFileName).c_str());}

	while(getline(modeFile, line)) {
		if(!line.empty()){
			tokens=line2Strings(line);
			if(idxPos < numAtoms && changeMode==false){isData=true;}
			else if(idxPos == numAtoms){isData=false;idxPos=0;}

			if(changeMode == true && tokens.size() > 1){
				modeIdx++;
				if (modeIdx == stoi(tokens.at(0)) && (modeIdx <= numModes)) {
					modeVal[modeIdx]=stof(tokens.at(1))*stof(tokens.at(1));
					eigVal[modeIdx-1]=0;
					//cout <<"# "<<modeVal[modeIdx]<< endl;
					changeMode=false;
				}
			}
			if(isData==true && tokens.size() == 3 ){
				float x=stof(tokens.at(0));
				float y=stof(tokens.at(1));
				float z=stof(tokens.at(2));

				eigVal[modeIdx-1]+=x*x+y*y+z*z;
				protModes[numModes*idxPos+modeIdx-1]						=stof(tokens.at(0));///sqrt(eigVal[modeIdx-1]);
				protModes[numAtoms*numModes+numModes*idxPos+modeIdx-1]		=stof(tokens.at(1));///sqrt(eigVal[modeIdx-1]);
				protModes[2*numAtoms*numModes+numModes*idxPos+modeIdx-1]	=stof(tokens.at(2));///sqrt(eigVal[modeIdx-1]);

				idxPos++;
				if(idxPos==numAtoms){changeMode=true;}
			}
		}
		idxLine++;
	}
	if (modeFile.bad()){	perror(("error while reading file " + modeFileName).c_str());}
	modeFile.close();

	modeX.insert(modeX.end(), modeY.begin(), modeY.end());
	modeX.insert(modeX.end(), modeZ.begin(), modeZ.end());
	for(int i=0;i<numAtoms;i++){
		for(int mode=0;mode<numModes;mode++){
			protModes[numModes*i+mode]						/=sqrt(eigVal[mode]);
			protModes[numAtoms*numModes+numModes*i+mode]		/=sqrt(eigVal[mode]);
			protModes[2*numAtoms*numModes+numModes*i+mode]	/=sqrt(eigVal[mode]);
		}
	}
	REAL* protBufMode = prot->getOrCreateModePtr();
	std::copy(protModes, protModes+ 3*numModes*numAtoms, protBufMode);
}
template<typename REAL>
void readGridFromGridFile(std::shared_ptr<GridUnion<REAL>> gridUnion, std::string filename) {
	/* original attract source used */
	attract::Grid grid;
	grid.read_std(filename.c_str());

	/* Inner grid */
	typename IntrplGrid<REAL>::Desc gridDescInner;
	memcpy(gridDescInner.typemask, grid.alphabet, 99 * sizeof(bool));
	gridDescInner.numGrids = grid.alphabetsize + 1;
	gridDescInner.dimN.x = grid.gridx;
	gridDescInner.dimN.y = grid.gridy;
	gridDescInner.dimN.z = grid.gridz;
	gridDescInner.gridSpacing = grid.gridspacing;
	gridDescInner.posMin.x = grid.ori[0];
	gridDescInner.posMin.y = grid.ori[1];
	gridDescInner.posMin.z = grid.ori[2];

	unsigned gridsize = gridDescInner.dimN.x * gridDescInner.dimN.y * gridDescInner.dimN.z
			* gridDescInner.numGrids;
	gridDescInner.grid = new float4[gridsize];
	memset(gridDescInner.grid, 0, gridsize * sizeof(float4));
	int gridpos = 0;
	for (int atomtype = 0; atomtype <= MAXATOMTYPES; atomtype++) {
		if (atomtype < MAXATOMTYPES && !grid.alphabet[atomtype])
			continue;
		for (unsigned int z = 0; z < gridDescInner.dimN.z; z++) {
			for (unsigned int y = 0; y < gridDescInner.dimN.y; y++) {
				for (unsigned int x = 0; x < gridDescInner.dimN.x; x++) {
					unsigned int index = gridDescInner.dimN.x * gridDescInner.dimN.y * z
							+ gridDescInner.dimN.x * y + x;
					attract::Potential &p = grid.innergrid[index].potential;
					int energradindex;
					energradindex = p[atomtype];
					if (energradindex > 0) {
						attract::EnerGradStd &e = grid.energrads_std[energradindex - 1];
						float4 &gg = gridDescInner.grid[gridpos];
						gg.x = e.grad[0];
						gg.y = e.grad[1];
						gg.z = e.grad[2];
						gg.w = e.energy;
					}
					gridpos++;
				}
			}
		}
	}

	/* put charge grid on first place. shift other grids by one */
	gridsize = gridDescInner.dimN.x * gridDescInner.dimN.y * gridDescInner.dimN.z;
	// copy charge grid
	float4* charge_grid = new float4[gridsize];
	std::copy(gridDescInner.grid+(gridDescInner.numGrids-1)*gridsize,
			gridDescInner.grid+gridDescInner.numGrids*gridsize,charge_grid);
	// shift other grids
	std::copy_backward(gridDescInner.grid, gridDescInner.grid+(gridDescInner.numGrids-1)*gridsize,
			gridDescInner.grid+gridDescInner.numGrids*gridsize);
	// copy charge grid to the beginning
	std::copy(charge_grid, charge_grid + gridsize, gridDescInner.grid);
	delete[] charge_grid;

	/* outer grid */
	typename IntrplGrid<REAL>::Desc gridDescOuter;
	memcpy(gridDescOuter.typemask, grid.alphabet, 99 * sizeof(bool));
	gridDescOuter.numGrids = grid.alphabetsize + 1;
	gridDescOuter.dimN.x = grid.gridx2;
	gridDescOuter.dimN.y = grid.gridy2;
	gridDescOuter.dimN.z = grid.gridz2;
	gridDescOuter.gridSpacing = 2 * grid.gridspacing;
	gridDescOuter.posMin.x = grid.ori[0] - grid.gridextension * grid.gridspacing;
	gridDescOuter.posMin.y = grid.ori[1] - grid.gridextension * grid.gridspacing;
	gridDescOuter.posMin.z = grid.ori[2] - grid.gridextension * grid.gridspacing;
	gridsize = gridDescOuter.dimN.x * gridDescOuter.dimN.y * gridDescOuter.dimN.z
			* gridDescOuter.numGrids;
	gridDescOuter.grid = new float4[gridsize];
	memset(gridDescOuter.grid, 0, gridsize * sizeof(float4));
	gridpos = 0;
	for (int atomtype = 0; atomtype <= MAXATOMTYPES; atomtype++) {
		if (atomtype < MAXATOMTYPES && !grid.alphabet[atomtype])
			continue;
		for (unsigned int z = 0; z < gridDescOuter.dimN.z; z++) {
			for (unsigned int y = 0; y < gridDescOuter.dimN.y; y++) {
				for (unsigned int x = 0; x < gridDescOuter.dimN.x; x++) {
					unsigned int index = gridDescOuter.dimN.x * gridDescOuter.dimN.y * z
							+ gridDescOuter.dimN.x * y + x;
					attract::Potential &p = grid.biggrid[index];
					int energradindex;
					energradindex = p[atomtype];
					if (energradindex > 0) {
						attract::EnerGradStd &e = grid.energrads_std[energradindex - 1];
						float4 &gg = gridDescOuter.grid[gridpos];
						gg.x = e.grad[0];
						gg.y = e.grad[1];
						gg.z = e.grad[2];
						gg.w = e.energy;
					}
					gridpos++;
				}
			}
		}
	}

	/* put charge grid on first place. shift other grids by one */
	gridsize = gridDescOuter.dimN.x * gridDescOuter.dimN.y * gridDescOuter.dimN.z;
	// copy charge grid
	charge_grid = new float4[gridsize];
	std::copy(gridDescOuter.grid+(gridDescOuter.numGrids-1)*gridsize,
			gridDescOuter.grid+gridDescOuter.numGrids*gridsize,charge_grid);
	// shift other grids
	std::copy_backward(gridDescOuter.grid, gridDescOuter.grid+(gridDescOuter.numGrids-1)*gridsize,
			gridDescOuter.grid+gridDescOuter.numGrids*gridsize);
	// copy charge grid to the beginning
	std::copy(charge_grid, charge_grid + gridsize, gridDescOuter.grid);
	delete[] charge_grid;

	/* NL grid */

	typename NLGrid<REAL>::Desc gridDescNL;
	gridDescNL.dimN.x = grid.gridx;
	gridDescNL.dimN.y = grid.gridy;
	gridDescNL.dimN.z = grid.gridz;
	gridDescNL.gridSpacing = grid.gridspacing;
	gridDescNL.dPlateau = grid.plateaudis;
	gridDescNL.posMin.x = grid.ori[0];
	gridDescNL.posMin.y = grid.ori[1];
	gridDescNL.posMin.z = grid.ori[2];
	gridDescNL.numEl = grid.nr_neighbours;
	gridDescNL.neighborArray = new uint[grid.nr_neighbours];
	for (int n = 0; n < grid.nr_neighbours; n++) {
		gridDescNL.neighborArray[n] = grid.neighbours[n].index;
	}
	gridsize = gridDescNL.dimN.x * gridDescNL.dimN.y * gridDescNL.dimN.z;
	gridDescNL.grid = new NeighbourDesc[gridsize];
	memset(gridDescNL.grid, 0, gridsize * sizeof(NeighbourDesc));
	for (unsigned int z = 0; z < gridDescNL.dimN.z; z++) {
		for (unsigned int y = 0; y < gridDescNL.dimN.y; y++) {
			for (unsigned int x = 0; x < gridDescNL.dimN.x; x++) {
				unsigned int index = gridDescNL.dimN.x * gridDescNL.dimN.y * z
						+ gridDescNL.dimN.x * y + x;
				attract::Voxel &v = grid.innergrid[index];
				NeighbourDesc &nd = gridDescNL.grid[index];
				nd.numEl = v.nr_neighbours;
				nd.idx = v.neighbourlist - 1; //in the ATTRACT code, 1 is the first neighbor index, not 0
			}
		}
	}
	// delete orig_grid resources
	delete [] grid.energrads_std;
	delete [] grid.energrads_torque;
	delete [] grid.neighbours;
	delete [] grid.biggrid;
	delete [] grid.innergrid;

	/* Create the GridUnion object */
	auto innerGrid = std::make_shared<IntrplGrid<REAL>>(gridDescInner);
	auto outerGrid = std::make_shared<IntrplGrid<REAL>>(gridDescOuter);
	auto NLgrid = std::make_shared<NLGrid<REAL>>(gridDescNL);

	gridUnion->inner = innerGrid;
	gridUnion->outer = outerGrid;
	gridUnion->NL = NLgrid;
	gridUnion->setTag(filename);
}

/*
 ** @brief: reads an ATTRACT forcefield parameter file and creates an
 ** ATTRACT parameter table object.
 */
template<typename REAL>
std::shared_ptr<ParamTable<REAL>> createParamTableFromFile(std::string filename) {
	// Check if REAL is of floating-point type
	auto table = std::make_shared<ParamTable<REAL>>();
	readParamTableFromFile(table, filename);
	return table;
}

/*
 ** @brief: as above but user provides AttrParamTable pointer.
 */
template<typename REAL>
void readParamTableFromFile(std::shared_ptr<ParamTable<REAL>> table, std::string filename) {
	using namespace std;

	ifstream infile(filename, ios::in);

	if (!infile.fail()) {
		unsigned potshape, numTypes;
		REAL swiOn, swiOff;
		infile >> potshape >> numTypes >> swiOn >> swiOff;


		table->setNumTypes(numTypes);
		table->setSwiOn(swiOn);
		table->setSwiOff(swiOff);

		typename ParamTable<REAL>::type_t* paramTable = table->getOrCreateTable();

		for (unsigned i = 0; i < numTypes; ++i) {
			for (unsigned j = 0; j < numTypes; ++j) {
				typename ParamTable<REAL>::type_t &params =paramTable[numTypes * i + j];
				infile >> params.rc;
			}
		}
		for (unsigned i = 0; i < numTypes; ++i) {
			for (unsigned j = 0; j < numTypes; ++j) {
				typename ParamTable<REAL>::type_t &params =paramTable[numTypes * i + j];
				infile >> params.ac;
			}
		}
		for (unsigned i = 0; i < numTypes; ++i) {
			for (unsigned j = 0; j < numTypes; ++j) {
				typename ParamTable<REAL>::type_t &params =paramTable[numTypes * i + j];
				infile >> params.ipon;
			}
		}

		infile.close();
		/* calculate rc, ac and other parameters */

		if (potshape == 8) {
			table->setPotShape(PotShape::_8_6);
			for (unsigned i = 0; i < numTypes; ++i) {
				for (unsigned j = 0; j < numTypes; ++j) {
					typename ParamTable<REAL>::type_t &params =paramTable[numTypes * i + j];
					REAL rbc = params.rc;
					REAL abc = params.ac;
					params.rc = abc * std::pow(rbc, potshape);

					params.ac = abc * std::pow(rbc, unsigned(6));
					if (params.ac > 0 && params.rc > 0) {
						params.emin = -27.0f * std::pow(params.ac, unsigned(4))
								/ (256.0f * std::pow(params.rc, unsigned(3)));
						params.rmin2 = 4.0f * params.rc
								/ (3.0f * params.ac);
					} else {
						params.emin = 0;
						params.rmin2 = 0;
					}
				}
			}
//
		} else if (potshape == 12) {
			table->setPotShape(PotShape::_12_6);
			for (unsigned i = 0; i < numTypes; ++i) {
				for (unsigned j = 0; j < numTypes; ++j) {
					typename ParamTable<REAL>::type_t &params =paramTable[numTypes * i + j];
					REAL rbc = params.rc;
					REAL abc = params.ac;
					params.rc = abc * pow(rbc, potshape);
					params.ac = abc * pow(rbc, unsigned(6));
					params.emin = -0.25 * abc;
					params.rmin2 = 1.25992105f * rbc * rbc;
				}
			}
		} else {
			infile.close();
			std::cerr << "Unknown potential shape " << potshape
					<< " in file: "<< filename << std::endl;
			exit(EXIT_FAILURE);
		}

	} else {
		infile.close();
		cerr << "Error: Failed to open file: " << filename << endl;
		exit(EXIT_FAILURE);
	}

}

std::vector<AttractEnGrad> readEnGradFromFile(std::string filename) {
//void readEnGradFromFile(std::string filename, unsigned& numEl) {
	using namespace std;

	ifstream file(filename);
	vector<AttractEnGrad> enGrads;

	string line;
	if (file.is_open()) {
		while (!file.eof()) {

			getline(file, line);
			if (!line.compare(0,2, "##")) {
				continue;
			} else if (!line.compare(0,8, " Energy:")){
				AttractEnGrad enGrad;
				{
					string substr = line.substr(8);
					stringstream stream(substr);
					stream >> enGrad.E;
				}

				{
					getline(file, line);
					stringstream stream(line);
					stream >> enGrad.E_VdW >> enGrad.E_El;
				}

				{
					getline(file, line);
					string substr = line.substr(11);
					stringstream stream(substr);
					stream >> enGrad.ang.x >> enGrad.ang.y >> enGrad.ang.z
						>> enGrad.pos.x >> enGrad.pos.y >> enGrad.pos.z;
				}
				enGrads.push_back(enGrad);
			}


		}
	} else {
		cerr << "Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	file.close();

	return enGrads;
}

std::vector<std::vector<DOF>> readDOF(std::string filename) {
	using namespace std;
	using namespace as;

	std::vector<std::vector<DOF>> DOF_molecules;

	ifstream file(filename);

	string line;
	size_t lineNo = 0;
	int i_molecules = 0;
	if (file.is_open()) {
		while (!file.eof()) {

			getline(file, line);
			++lineNo;


			if (!line.compare(0,1, "#")) { // 0 == true
				continue;
			}

			/* read all dofs until the next "#" */
			unsigned i = 0;
			while (line.compare(0,1, "#") != 0 && !file.eof()) {

				if (i_molecules == 0) {
					DOF_molecules.push_back(std::vector<DOF> ());
				}

				std::vector<DOF>& vec = DOF_molecules[i];
				DOF dof ;
				{
					stringstream stream(line);
					stream >> dof._6D.ang.x >> dof._6D.ang.y >> dof._6D.ang.z
						>> dof._6D.pos.x >> dof._6D.pos.y >> dof._6D.pos.z;
					for (size_t i = 0; !stream.eof(); ++i) {
						stream >> dof.dofs[i];
						if(stream.fail()) {
							stream.clear();
							string parsed;
							stream >> parsed;
							errorDOFFormat(filename);
							cerr << "Cannot parse '" << parsed << "' at #" << i_molecules << " (line " << lineNo << ")" << endl;
							exit(EXIT_FAILURE);
						}
					}
				}
				vec.push_back(dof);

				++i;
				getline(file, line);
				++lineNo;
			}
			/* check if i equals the number of molecules == DOF_molecules.size(),
			 * otherwise we miss a molecule in the definition */
			if (i != DOF_molecules.size()) {
				errorDOFFormat(filename);
				cerr << "The DOF definition is incomplete at #" << i_molecules << " (line " << lineNo << ")" << endl;
						exit(EXIT_FAILURE);
			}
			++i_molecules;
		}
	} else {
		cerr << "Error: Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	file.close();

	return DOF_molecules;

}

template<typename REAL>
DOFHeader<REAL> readDOFHeader(std::string filename) {
	using namespace std;
	using vec3_t = Vec3<REAL>;
	ifstream file(filename);
	string line;

	DOFHeader<REAL> header;

	bool finished_reading_pivots = false;
	bool finished_reading_centered_receptor = false;
	bool finished_reading_centered_ligands = false;
	header.auto_pivot = false;

	if (file.is_open()) {
		while (!file.eof()) {

			if (finished_reading_pivots && finished_reading_centered_receptor && finished_reading_centered_ligands) {
				break;
			}

			getline(file, line);

			/* Check for commandline comment (use in ATTRACT output) */
			if (!line.compare(0,2, "##")) { // 0 == true
				continue;
			}

			/* get pivots */

			/* split line to vector of strings */
			vector<string> tokens = line2Strings(line);
//			cout << line << endl;
			if (tokens[0].compare("#pivot") == 0 && !finished_reading_pivots) {
				if (tokens[1].compare("auto") == 0) {
					header.auto_pivot = true;
					finished_reading_pivots = true;
					continue;
				} else {
					if (tokens.size() != 5) { // validate format
						errorDOFFormat(filename);
						cerr << "Incorrect pivot definition." << endl;
						exit(EXIT_FAILURE);
					} else { // read pivots
						vec3_t pivot;
						stringstream(tokens[2]) >> pivot.x;
						stringstream(tokens[3]) >> pivot.y;
						stringstream(tokens[4]) >> pivot.z;
						header.pivots.push_back(pivot);
						continue;
					}
				}

			} else if (tokens[0].compare("#pivot") == 0 && finished_reading_pivots) { // validate format
				errorDOFFormat(filename);
				cerr << "Incorrect pivot definition." << endl;
				exit(EXIT_FAILURE);
			}

			/* read centered definition */
			if (tokens[0].compare("#centered") == 0 &&
					(!finished_reading_centered_receptor || !finished_reading_centered_ligands)) {
				finished_reading_pivots = true;
				if (tokens[1].compare("receptor:") == 0) {
					if (tokens[2].compare("true") == 0) {
						header.centered_receptor = true;
						finished_reading_centered_receptor = true;
						continue;
					} else if (tokens[2].compare("false") == 0) {
						header.centered_receptor = false;
						finished_reading_centered_receptor = true;
						continue;
					} else {
						errorDOFFormat(filename);
						cerr << "Receptor centering definition invalid." << endl;
						exit(EXIT_FAILURE);
					}
				} else if (tokens[1].compare("ligand:") == 0 || tokens[1].compare("ligands:") == 0) {
					if (tokens[2].compare("true") == 0) {
						header.centered_ligands = true;
						finished_reading_centered_ligands = true;
						continue;
					} else if (tokens[2].compare("false") == 0) {
						header.centered_ligands = false;
						finished_reading_centered_ligands = true;
						continue;
					} else {
						errorDOFFormat(filename);
						cerr << "Ligand centering definition invalid." << endl;
						exit(EXIT_FAILURE);
					}
				} else {
					errorDOFFormat(filename);
					cerr << "Centering definition invalid." << endl;
					exit(EXIT_FAILURE);
				}
			}

		} // while (!file.eof())
	} else {
		cerr << "Error: Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	return header;
}

std::vector<std::string> readFileNamesFromEnsembleList(std::string filename) {
	using namespace std;
	using namespace as;

	ifstream file(filename);

	vector<string> fileNames;
	string line;
	if (file.is_open()) {
		getline(file, line);
		while (!file.eof()) {
			fileNames.push_back(line);
			getline(file, line);
		}
	} else {
		cerr << "Error: Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	file.close();
	return fileNames;
}

static void errorGridAlphabetFormat(std::string filename) {
	std::cerr << "Error reading file " << filename
			<< ". Grid alphabet format is not supported." << std::endl;
}

std::vector<unsigned> readGridAlphabetFromFile(std::string filename) {

	using namespace std;
	using namespace as;

	ifstream file(filename);
	vector<unsigned> vec;
	string line;
	if (file.is_open()) {
		while (std::getline(file, line))
		{
			std::istringstream iss(line);
			int key = 0; // serves as a key value pair for a map
			if (!(iss >> key) || key < 0) {
				errorGridAlphabetFormat(filename);
				exit(EXIT_FAILURE);
			}
			vec.push_back(static_cast<unsigned>(key));
		}
	} else {
		cerr << "Error: Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	return vec;
}


template
std::shared_ptr<Protein<float>> createProteinFromPDB(std::string filename);

template
void readProteinFromPDB(std::shared_ptr<Protein<float>> prot, std::string filename);

template
void readHMMode(std::shared_ptr<Protein<float>> prot, std::string modeFileName);

template
void readHMMode(std::shared_ptr<Protein<double>> prot, std::string modeFileName);

template
std::shared_ptr<Protein<double>> createProteinFromPDB(std::string filename);

template
void readProteinFromPDB(std::shared_ptr<Protein<double>> prot, std::string filename);

template
std::shared_ptr<GridUnion<float>> createGridFromGridFile(std::string filename);

template
void readGridFromGridFile(std::shared_ptr<GridUnion<float>>, std::string filename);

template
std::shared_ptr<GridUnion<double>> createGridFromGridFile(std::string filename);

template
void readGridFromGridFile(std::shared_ptr<GridUnion<double>>, std::string filename);

template
std::shared_ptr<ParamTable<float>> createParamTableFromFile(std::string filename);

template
void readParamTableFromFile(std::shared_ptr<ParamTable<float>>, std::string filename);

template
std::shared_ptr<ParamTable<double>> createParamTableFromFile(std::string filename);

template
void readParamTableFromFile(std::shared_ptr<ParamTable<double>>, std::string filename);

template
DOFHeader<float> readDOFHeader<float>(std::string filename);

template
DOFHeader<double> readDOFHeader<double>(std::string filename);

} // namespace as



//unsigned asDB::readProteinSizeFromPDB(std::string filename) {
//	using namespace std;
//
//	ifstream file(filename);
//	string line;
//
//	unsigned count = 0;
//
//	if (file.is_open()) {
//		getline(file, line);
//		while (!file.eof()) {
//
//			/* split line to vector of strings */
//			vector<string> tokens = line2Strings(line);
//
//			if (!(tokens.size() == 12 || tokens.size() == 13)
//					&& !line.empty()) { // 12: old reduced format, 13: new reduced format
//				errorPDBFormat(filename);
//				cerr << "Expected 12 or 13 columns." << endl;
//				exit(EXIT_FAILURE);
//			}
//			/* empty lines are acceptable, e.g. at the end */
//			if (line.compare(0, 4, "ATOM") != 0 && !line.empty()) { // does not match (0==equal)
//				errorPDBFormat(filename);
//				cerr << "Expected 'ATOM' to be in the first column." << endl;
//				exit(EXIT_FAILURE);
//			}
//
//			if (line.compare(0, 4, "ATOM") == 0) {
//				++count;
//			}
//
//			getline(file, line);
//		}
//	} else {
//		cerr << "Failed to open file " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//	file.close();
//
//	return count;
//}
//
//Protein* asDB::createProteinFromDumpFile(std::string filename) {
//	using namespace std;
//
////	string path = "../data/";
////	string path = "/home/uwe/nsight/cuda-workspace/AttractServer/data/";
////	string filename = path + name;
////	cout << filename << endl;
////	string filename = name;
//
//	ifstream file(filename.c_str(), ios::in | ios::binary);
//	if(!file.is_open()) {
//		cerr << "Error: Failed to open file " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	ProteinDesc desc;
//
//	if (!file.read((char*) &desc.numAtoms, sizeof(unsigned))) {
//		cerr << "Error read: numAtoms" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	float *pos = new float[3 * desc.numAtoms];
//	unsigned *type = new unsigned[desc.numAtoms];
//	float *charge = new float[desc.numAtoms];
//
//	if (!file.read((char*) pos, 3 * desc.numAtoms * sizeof(float))) {
//		cerr << "Error read: pos" << endl;
//		exit(EXIT_FAILURE);
//	}
//	desc.pos = pos;
//
//	if (!file.read((char*) type, desc.numAtoms * sizeof(unsigned))) {
//		cerr << "Error read: type" << endl;
//		exit(EXIT_FAILURE);
//	}
//	desc.type = type;
//
//	if (!file.read((char*) charge, desc.numAtoms * sizeof(float))) {
//		cerr << "Error read: charge" << endl;
//		exit(EXIT_FAILURE);
//	}
//	desc.charge = charge;
//
//	if (!file.read((char*) &desc.numModes, sizeof(unsigned))) {
//		cerr << "Error read: numModes" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (desc.numModes > 0) {
//		float *modes = new float[3 * desc.numAtoms * desc.numModes];
//		if (!file.read((char*) modes,
//				3 * desc.numAtoms * desc.numModes * sizeof(float))) {
//			cerr << "Error read: modes" << endl;
//			exit(EXIT_FAILURE);
//		}
//		desc.modes = modes;
//	} else {
//		desc.modes = nullptr;
//	}
//
//	file.close();
//
//	Protein* protein = new Protein(desc);
//	protein->setTag(filename);
//
//	return protein;
//}
//
//GridUnion* asDB::createGridFromGridFile(std::string filename) {
//	GridUnion* grid = new GridUnion();
//	readGridFromGridFile(grid, filename);
//	return grid;
//}
//
//void asDB::readGridFromGridFile(GridUnion* gridUnion, std::string filename) {
//
//	/* original attract source used */
//	Grid grid;
//	grid.read_std(filename.c_str());
//
//	/* Inner grid */
//	GradEnGridDesc gridDescInner;
//	memcpy(gridDescInner.typemask, grid.alphabet, 99 * sizeof(bool));
//	gridDescInner.numGrids = grid.alphabetsize + 1;
//	gridDescInner.width = grid.gridx;
//	gridDescInner.height = grid.gridy;
//	gridDescInner.depth = grid.gridz;
//	gridDescInner.gridSpacing = grid.gridspacing;
//	for (int n = 0; n < 3; n++)
//		gridDescInner.posMin[n] = grid.ori[n];
//	unsigned int gridsize = gridDescInner.width * gridDescInner.height * gridDescInner.depth
//			* gridDescInner.numGrids;
//	gridDescInner.grid = new float4[gridsize];
//	memset(gridDescInner.grid, 0, gridsize * sizeof(float4));
//	int gridpos = 0;
//	for (int atomtype = 0; atomtype <= MAXATOMTYPES; atomtype++) {
//		if (atomtype < MAXATOMTYPES && !grid.alphabet[atomtype])
//			continue;
//		for (unsigned int z = 0; z < gridDescInner.depth; z++) {
//			for (unsigned int y = 0; y < gridDescInner.height; y++) {
//				for (unsigned int x = 0; x < gridDescInner.width; x++) {
//					unsigned int index = gridDescInner.width * gridDescInner.height * z
//							+ gridDescInner.width * y + x;
//					Potential &p = grid.innergrid[index].potential;
//					int energradindex;
//					energradindex = p[atomtype];
//					if (energradindex > 0) {
//						EnerGradStd &e = grid.energrads_std[energradindex - 1];
//						float4 &gg = gridDescInner.grid[gridpos];
//						gg.x = e.grad[0];
//						gg.y = e.grad[1];
//						gg.z = e.grad[2];
//						gg.w = e.energy;
//					}
//					gridpos++;
//				}
//			}
//		}
//	}
//
////	print(gridDescInner, 31, 31 , 0,0, 25, 30, 25, 30, 25, 30);
//
//	/* put charge grid on first place. shift other grids by one */
//	if (1) {
//		unsigned gridsize = gridDescInner.width * gridDescInner.height * gridDescInner.depth;
//		// copy charge grid
//		float4* charge_grid = new float4[gridsize];
//		std::copy(gridDescInner.grid+(gridDescInner.numGrids-1)*gridsize,
//				gridDescInner.grid+gridDescInner.numGrids*gridsize,charge_grid);
//		// shift other grids
//		std::copy_backward(gridDescInner.grid, gridDescInner.grid+(gridDescInner.numGrids-1)*gridsize,
//				gridDescInner.grid+gridDescInner.numGrids*gridsize);
//		// copy charge grid to the beginning
//		std::copy(charge_grid, charge_grid + gridsize, gridDescInner.grid);
//		delete[] charge_grid;
//	}
//
//	/* outer grid */
//	GradEnGridDesc gridDescOuter;
//	memcpy(gridDescOuter.typemask, grid.alphabet, 99 * sizeof(bool));
//	gridDescOuter.numGrids = grid.alphabetsize + 1;
//	gridDescOuter.width = grid.gridx2;
//	gridDescOuter.height = grid.gridy2;
//	gridDescOuter.depth = grid.gridz2;
//	gridDescOuter.gridSpacing = 2 * grid.gridspacing;
//	for (int n = 0; n < 3; n++)
//		gridDescOuter.posMin[n] = grid.ori[n]
//				- grid.gridextension * grid.gridspacing;
//	gridsize = gridDescOuter.width * gridDescOuter.height * gridDescOuter.depth
//			* gridDescOuter.numGrids;
//	gridDescOuter.grid = new float4[gridsize];
//	memset(gridDescOuter.grid, 0, gridsize * sizeof(float4));
//	gridpos = 0;
//	for (int atomtype = 0; atomtype <= MAXATOMTYPES; atomtype++) {
//		if (atomtype < MAXATOMTYPES && !grid.alphabet[atomtype])
//			continue;
//		for (unsigned int z = 0; z < gridDescOuter.depth; z++) {
//			for (unsigned int y = 0; y < gridDescOuter.height; y++) {
//				for (unsigned int x = 0; x < gridDescOuter.width; x++) {
//					unsigned int index = gridDescOuter.width * gridDescOuter.height * z
//							+ gridDescOuter.width * y + x;
//					Potential &p = grid.biggrid[index];
//					int energradindex;
//					energradindex = p[atomtype];
//					if (energradindex > 0) {
//						EnerGradStd &e = grid.energrads_std[energradindex - 1];
//						float4 &gg = gridDescOuter.grid[gridpos];
//						gg.x = e.grad[0];
//						gg.y = e.grad[1];
//						gg.z = e.grad[2];
//						gg.w = e.energy;
//					}
//					gridpos++;
//				}
//			}
//		}
//	}
//
//	/* put charge grid on first place. shift other grids by one */
//	if (1) {
//		unsigned gridsize = gridDescOuter.width * gridDescOuter.height * gridDescOuter.depth;
//		// copy charge grid
//		float4* charge_grid = new float4[gridsize];
//		std::copy(gridDescOuter.grid+(gridDescOuter.numGrids-1)*gridsize,
//				gridDescOuter.grid+gridDescOuter.numGrids*gridsize,charge_grid);
//		// shift other grids
//		std::copy_backward(gridDescOuter.grid, gridDescOuter.grid+(gridDescOuter.numGrids-1)*gridsize,
//				gridDescOuter.grid+gridDescOuter.numGrids*gridsize);
//		// copy charge grid to the beginning
//		std::copy(charge_grid, charge_grid + gridsize, gridDescOuter.grid);
//		delete[] charge_grid;
//	}
//
//	/* NL grid */
//
//	NLGridDesc gridDescNL;
//	gridDescNL.width = grid.gridx;
//	gridDescNL.height = grid.gridy;
//	gridDescNL.depth = grid.gridz;
//	gridDescNL.gridSpacing = grid.gridspacing;
//	gridDescNL.dPlateau = grid.plateaudis;
//	for (int n = 0; n < 3; n++)
//		gridDescNL.posMin[n] = grid.ori[n];
//	gridDescNL.numEl = grid.nr_neighbours;
//	gridDescNL.neighborArray = new uint[grid.nr_neighbours];
//	for (int n = 0; n < grid.nr_neighbours; n++) {
//		gridDescNL.neighborArray[n] = grid.neighbours[n].index;
//	}
//	gridsize = gridDescNL.width * gridDescNL.height * gridDescNL.depth;
//	gridDescNL.grid = new NeighbourDesc[gridsize];
//	memset(gridDescNL.grid, 0, gridsize * sizeof(NeighbourDesc));
//	for (unsigned int z = 0; z < gridDescNL.depth; z++) {
//		for (unsigned int y = 0; y < gridDescNL.height; y++) {
//			for (unsigned int x = 0; x < gridDescNL.width; x++) {
//				unsigned int index = gridDescNL.width * gridDescNL.height * z
//						+ gridDescNL.width * y + x;
//				Voxel &v = grid.innergrid[index];
//				NeighbourDesc &nd = gridDescNL.grid[index];
//				nd.numEl = v.nr_neighbours;
//				nd.idx = v.neighbourlist - 1; //in the ATTRACT code, 1 is the first neighbor index, not 0
//			}
//		}
//	}
//
//	/* Create the GridUnion object */
//	IntrplGrid* innerGrid = new IntrplGrid(gridDescInner);
//	IntrplGrid* outerGrid = new IntrplGrid(gridDescOuter);
//	NLGrid* NLgrid = new NLGrid(gridDescNL);
//
//	gridUnion->setInnerGrid(innerGrid);
//	gridUnion->setOuterGrid(outerGrid);
//	gridUnion->setNLGrid(NLgrid);
//	gridUnion->setTag(filename);
//
//}
//
//GridUnion* asDB::createGridUnionFromDumpFile(std::string filename) {
//	using namespace std;
//	//	string path = "../data/";
////	string path = "/home/uwe/nsight/cuda-workspace/AttractServer/data/";
////	string filename = path + name;
////	cout << filename << endl;
////	string filename = name;
//	ifstream file(filename.c_str(), ios::in | ios::binary);
//	if(!file.is_open()) {
//		cerr << "Error: Failed to open file " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	/* reading the gradient-energy grids */
//	GradEnGridDesc intrplDesc[2];
//	for (int k = 0; k < 2; ++k) {
//		if (!file.read((char*) intrplDesc[k].typemask, 99 * sizeof(bool))) {
//			cerr << "Error read: typemask" << endl;
//			exit(EXIT_FAILURE);
//		}
//
//		/* Analyze typemask */
//		unsigned numGrids = 0;
//		for (unsigned i = 0; i < 99; ++i) {
//			if (intrplDesc[k].typemask[i] == true) {
//				++numGrids;
//			}
//		}
//		/* Number is incremented by one to count for the charge Grid */
//		intrplDesc[k].numGrids = numGrids + 1;
//
//		if (!file.read((char*) &intrplDesc[k].width, sizeof(unsigned))) {
//			cerr << "Error read: width" << endl;
//			exit(EXIT_FAILURE);
//		}
//		if (!file.read((char*) &intrplDesc[k].height, sizeof(unsigned))) {
//			cerr << "Error read: height" << endl;
//			exit(EXIT_FAILURE);
//		}
//		if (!file.read((char*) &intrplDesc[k].depth, sizeof(unsigned))) {
//			cerr << "Error read: depth" << endl;
//			exit(EXIT_FAILURE);
//		}
//		if (!file.read((char*) &intrplDesc[k].gridSpacing, sizeof(float))) {
//			cerr << "Error read: gridSpacing" << endl;
//			exit(EXIT_FAILURE);
//		}
//		if (!file.read((char*) intrplDesc[k].posMin, 3 * sizeof(float))) {
//			cerr << "Error read: posMin" << endl;
//			exit(EXIT_FAILURE);
//		}
//
//		unsigned gridsize = intrplDesc[k].width * intrplDesc[k].height * intrplDesc[k].depth;
//
//		unsigned numEl = intrplDesc[k].numGrids * gridsize;
//		float4* grid = new float4[numEl];
//
////		if (!file.read((char*) grid, numEl * sizeof(float4))) {
////			cerr << "Error read: grid" << endl;
////			exit(EXIT_FAILURE);
////		}
//
//		/* first, read not the charge grid -> leaf space for charge grid */
//		if (!file.read((char*) (grid + gridsize), (numEl-gridsize) * sizeof(float4))) {
//			cerr << "Error read: grid" << endl;
//			exit(EXIT_FAILURE);
//		}
//
//		/* second, read the charge grid */
//		if (!file.read((char*) grid, gridsize * sizeof(float4))) {
//			cerr << "Error read: grid" << endl;
//			exit(EXIT_FAILURE);
//		}
//
//		intrplDesc[k].grid = grid;
//	}
//
//	/* reading the neighbor list grid */
//	NLGridDesc NLdesc;
//
//	if (!file.read((char*) &NLdesc.width, sizeof(unsigned))) {
//		cerr << "Error read: width" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (!file.read((char*) &NLdesc.height, sizeof(unsigned))) {
//		cerr << "Error read: height" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (!file.read((char*) &NLdesc.depth, sizeof(unsigned))) {
//		cerr << "Error read: depth" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (!file.read((char*) &NLdesc.gridSpacing, sizeof(float))) {
//		cerr << "Error read: gridSpacing" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (!file.read((char*) &NLdesc.dPlateau, sizeof(float))) {
//		cerr << "Error read: dPlateau" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if (!file.read((char*) NLdesc.posMin, 3 * sizeof(float))) {
//		cerr << "Error read: posMin" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	unsigned numEl = NLdesc.width * NLdesc.height * NLdesc.depth;
//	NeighbourDesc* grid = new NeighbourDesc[numEl];
//
//	if (!file.read((char*) grid, numEl * sizeof(NeighbourDesc))) {
//		cerr << "Error read: grid" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	NLdesc.grid = grid;
//
//	if (!file.read((char*) &NLdesc.numEl, sizeof(unsigned))) {
//		cerr << "Error read: width" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	numEl = NLdesc.numEl;
//	unsigned* neighbourArray = new unsigned[numEl];
//
//	if (!file.read((char*) neighbourArray, numEl * sizeof(unsigned))) {
//		cerr << "Error read: grid" << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	NLdesc.neighborArray = neighbourArray;
//
//	/* Create the GridUnion object */
//	IntrplGrid* innerGrid = new IntrplGrid(intrplDesc[0]);
//	IntrplGrid* outerGrid = new IntrplGrid(intrplDesc[1]);
//	NLGrid* NLgrid = new NLGrid(NLdesc);
//
//	GridUnion* gridUnion = new GridUnion(innerGrid, outerGrid,
//			NLgrid);
//	gridUnion->setTag(filename);
//
//	file.close();
//
//	return gridUnion;
//}
//
//AttrParamTable* asDB::createParamTableFromFile(std::string filename) {
//	AttrParamTable* table = new AttrParamTable;
//	readParamTableFromFile(table, filename);
//	return table;
//}
//
//void asDB::readParamTableFromFile(AttrParamTable* table, std::string filename) {
//	using namespace std;
//	using namespace as;
//
//	ifstream infile(filename, ios::in);
//
//	if (!infile.fail()) {
//		unsigned potshape, numTypes;
//		float swiOn, swiOff;
//		infile >> potshape >> numTypes >> swiOn >> swiOff;
//
//
//		table->setNumTypes(numTypes);
//		table->setSwiOn(swiOn);
//		table->setSwiOff(swiOff);
//
//		AttrParamTable::type* paramTable = table->getOrCreateTable();
//
//		for (unsigned i = 0; i < numTypes; ++i) {
//			for (unsigned j = 0; j < numTypes; ++j) {
//				AttrParamTable::type &params =paramTable[numTypes * i + j];
//				infile >> params.rc;
//			}
//		}
//		for (unsigned i = 0; i < numTypes; ++i) {
//			for (unsigned j = 0; j < numTypes; ++j) {
//				AttrParamTable::type &params =paramTable[numTypes * i + j];
//				infile >> params.ac;
//			}
//		}
//		for (unsigned i = 0; i < numTypes; ++i) {
//			for (unsigned j = 0; j < numTypes; ++j) {
//				AttrParamTable::type &params =paramTable[numTypes * i + j];
//				infile >> params.ipon;
//			}
//		}
//
//		infile.close();
//		/* calculate rc, ac and other parameters */
//
//		if (potshape == 8) {
//			table->setPotShape(AttrParamTable::PotShape::_8_6);
//			for (unsigned i = 0; i < numTypes; ++i) {
//				for (unsigned j = 0; j < numTypes; ++j) {
//					AttrParamTable::type &params =paramTable[numTypes * i + j];
//					float rbc = params.rc;
//					float abc = params.ac;
//					params.rc = abc * std::pow(rbc, potshape);
//
//					params.ac = abc * std::pow(rbc, unsigned(6));
//					if (params.ac > 0 && params.rc > 0) {
//						params.emin = -27.0f * std::pow(params.ac, unsigned(4))
//								/ (256.0f * std::pow(params.rc, unsigned(3)));
//						params.rmin2 = 4.0f * params.rc
//								/ (3.0f * params.ac);
//					} else {
//						params.emin = 0;
//						params.rmin2 = 0;
//					}
//				}
//			}
////
//		} else if (potshape == 12) {
//			table->setPotShape(AttrParamTable::PotShape::_12_6);
//			for (unsigned i = 0; i < numTypes; ++i) {
//				for (unsigned j = 0; j < numTypes; ++j) {
//					AttrParamTable::type &params =paramTable[numTypes * i + j];
//					float rbc = params.rc;
//					float abc = params.ac;
//					params.rc = abc * pow(rbc, potshape);
//					params.ac = abc * pow(rbc, unsigned(6));
//					params.emin = -0.25 * abc;
//					params.rmin2 = 1.25992105f * rbc * rbc;
//				}
//			}
//		} else {
//			infile.close();
//			std::cerr << "Unknown potential shape " << potshape
//					<< " in file: "<< filename << std::endl;
//			exit(EXIT_FAILURE);
//		}
//
//	} else {
//		infile.close();
//		cerr << "Error: Failed to open file: " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//
//}
//
//
//void asDB::readEnsembleDOFFromFile(std::string filename, std::vector<std::vector<DOF>>& DOF_molecules ) {
//	using namespace std;
//	using namespace as;
//
//	ifstream file(filename);
//
//	string line;
//	int i_molecules = 0;
//	if (file.is_open()) {
//		while (!file.eof()) {
//
//			getline(file, line);
//
//
//			if (!line.compare(0,1, "#")) { // 0 == true
//				continue;
//			}
//
//			/* read all dofs until the next "#" */
//			unsigned i = 0;
//			while (line.compare(0,1, "#") != 0 && !file.eof()) {
//
//				if (i_molecules == 0) {
//					DOF_molecules.push_back(std::vector<DOF> ());
//				}
//
//				std::vector<DOF>& vec = DOF_molecules[i];
//				DOF dof ;
//				{
//					stringstream stream(line);
//					stream >> dof.ligId >> dof.ang.x >> dof.ang.y >> dof.ang.z
//						>> dof.pos.x >> dof.pos.y >> dof.pos.z;
//				}
//				vec.push_back(dof);
//
//				++i;
//				getline(file, line);
//			}
//			/* check if i equals the number of molecules == DOF_molecules.size(),
//			 * otherwise we miss a molecule in the definition */
//			if (i != DOF_molecules.size()) {
//				errorDOFFormat(filename);
//				cerr << "The DOF definition is incomplete at #" << i_molecules << "." << endl;
//						exit(EXIT_FAILURE);
//			}
//			++i_molecules;
//		}
//	} else {
//		cerr << "Error: Failed to open file " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	file.close();
//
//}
//

//
//std::vector<std::string> asDB::readFileNamesFromEnsembleList(std::string filename) {
//	using namespace std;
//	using namespace as;
//
//	ifstream file(filename);
//
//	vector<string> fileNames;
//	string line;
//	if (file.is_open()) {
//		getline(file, line);
//		while (!file.eof()) {
//			fileNames.push_back(line);
//			getline(file, line);
//		}
//	} else {
//		cerr << "Error: Failed to open file " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	file.close();
//	return fileNames;
//}
//
//static void errorGridAlphabetFormat(std::string filename) {
//	std::cerr << "Error reading file " << filename
//			<< ". Grid alphabet format is not supported." << std::endl;
//}
//
//std::vector<unsigned> asDB::readGridAlphabetFromFile(std::string filename) {
//
//	using namespace std;
//	using namespace as;
//
//	ifstream file(filename);
//	vector<unsigned> vec;
//	string line;
//	if (file.is_open()) {
//		while (std::getline(file, line))
//		{
//			std::istringstream iss(line);
//			int key = 0; // serves as a key value pair for a map
//			if (!(iss >> key) || key < 0) {
//				errorGridAlphabetFormat(filename);
//				exit(EXIT_FAILURE);
//			}
//			vec.push_back(static_cast<unsigned>(key));
//		}
//	} else {
//		cerr << "Error: Failed to open file " << filename << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	return vec;
//}

