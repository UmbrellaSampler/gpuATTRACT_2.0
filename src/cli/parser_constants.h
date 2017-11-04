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

/**
 * Two body docking context
 */

constexpr auto DEFAULT_RECEPTOR_PDB_FILE = "receptorr.pdb";
constexpr auto DEFAULT_LIGANG_PDB_FILE = "ligandr.pdb";
constexpr auto DEFAULT_RECEPTOR_GRID_FILE = "receptorgrid.grid";
constexpr auto DEFAULT_PARAMETER_FILE = "attract.par";
constexpr auto DEFAULT_GRID_ALPAHBET_FILE = "receptorgrid.alphabet";

/**
 * Server context
 */
constexpr int DEFAULT_NUM_CPUS = 1;

constexpr int DEFAULT_CHUNK_SIZE = 1000;

/**
 * Simulation context
 */
constexpr auto DEFAULT_DIELEC_MODE = "variable";
constexpr double DEFAULT_EPSILON_CONSTANT = 15.0;

#endif /* SRC_CLI_PARSER_CONSTANTS_H_ */
