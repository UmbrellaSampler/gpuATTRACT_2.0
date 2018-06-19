/*
 * parser_constants.h
 *
 *  Created on: Jul 23, 2017
 *      Author: uwe
 */

#ifndef SRC_CLI_PARSER_CONSTANTS_H_
#define SRC_CLI_PARSER_CONSTANTS_H_

#include <string>

/**
 * Application context
 */
constexpr auto APP_SHORT_NAME_SC = "sc";
constexpr auto APP_SHORT_NAME_EM = "em";

/**
 * Two body docking context
 */


/* Input Files */
constexpr auto FILE_DEFAULT_RECEPTOR_PDB = "receptorr.pdb";
constexpr auto FILE_DEFAULT_LIGANG_PDB = "ligandr.pdb";
constexpr auto FILE_DEFAULT_RECEPTOR_GRID = "receptorgrid.grid";
constexpr auto FILE_DEFAULT_LIGAND_GRID = "ligandgrid.grid";
constexpr auto FILE_DEFAULT_PARAMETER = "attract.par";
constexpr auto FILE_DEFAULT_GRID_ALPAHBET_RECEPTOR = "receptorgrid.alphabet";
constexpr auto FILE_DEFAULT_GRID_ALPAHBET_LIGAND = "ligandgrid.alphabet";
constexpr auto DEFAULT_MODE_RECEPTOR_FILE = "modesReceptor.dat";
constexpr auto DEFAULT_MODE_LIGAND_FILE = "modesLigand.dat";

/* Generic */
constexpr auto GENERIC_DEFAULT_PRECISION = "single";
constexpr auto GENERIC_ALLOWED_PRECISION = {"single", "double"};

/* Simulation */
constexpr auto SIM_DEFAULT_DIELEC = "variable";
constexpr auto SIM_ALLOWED_DIELEC = {"variable", "constant"};
constexpr auto SIM_DEFAULT_EPSILON = 15.0;
constexpr int  DEFAULT_NUM_MODES = 0;
constexpr double DEFAULT_CUTOFF = -1.0;
constexpr double DEFAULT_MODEFORCEFAC = 1.0;

/**
 * Server context
 */

constexpr int SERVER_DEFAULT_NUM_CPUS = 1;
constexpr int SERVER_DEFAULT_CHUNK_SIZE = 1000;
constexpr int SERVER_DEFAULT_DEVICE_ID = 0;

/**
 * EM
 */
constexpr auto EM_DEFAULT_SOLVER = "VA13";
constexpr auto EM_ALLOWED_SOLVERS = {"VA13", "BFGS"};
constexpr auto EM_DEFAULT_CONCURRENCY = 20000;
constexpr auto EM_DEFAULT_NUM_CHUNKS = 2;



/**
 * Simulation context
 */

#endif /* SRC_CLI_PARSER_CONSTANTS_H_ */
