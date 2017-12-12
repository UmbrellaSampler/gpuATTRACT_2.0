/*
 * DOFConverter.h
 *
 *  Created on: Dec 12, 2017
 *      Author: glenn
 */

#ifndef DOFCONVERTER_H_
#define DOFCONVERTER_H_
#include "Types_6D_Modes.h"
#include "Types_6D.h"
#include "readFile.h"
#include <vector>

namespace as{
template<typename REAL>
std::vector<std::vector<DOF_6D<REAL>>> DOFConverter(std::vector<std::vector<DOF>> &in) {
	std::vector<DOF_6D<REAL>> outRec(in[0].size());
	std::vector<DOF_6D<REAL>> outLig(in[1].size());

		for (int i=0; i < in[0].size(); i++){

			outRec[i].ang.x = in[0][i]._6D.ang.x;
			outRec[i].ang.y = in[0][i]._6D.ang.y;
			outRec[i].ang.z = in[0][i]._6D.ang.z;
			outRec[i].pos.x = in[0][i]._6D.pos.x;
			outRec[i].pos.y = in[0][i]._6D.pos.y;
			outRec[i].pos.z = in[0][i]._6D.pos.z;
		}

		for (int i=0; i < in[1].size(); i++){
			outLig[i].ang.x = in[1][i]._6D.ang.x;
			outLig[i].ang.y = in[1][i]._6D.ang.y;
			outLig[i].ang.z = in[1][i]._6D.ang.z;
			outLig[i].pos.x = in[1][i]._6D.pos.x;
			outLig[i].pos.y = in[1][i]._6D.pos.y;
			outLig[i].pos.z = in[1][i]._6D.pos.z;
		}
		std::vector<std::vector<DOF_6D<REAL>>> outVec;
		outVec.push_back(outRec);
		outVec.push_back(outLig);
		return outVec;
}


template<typename REAL>
std::vector<std::vector<DOF_6D_Modes<REAL>>> DOFConverter_Modes(std::vector<std::vector<DOF>>  &in) {
	std::vector<DOF_6D_Modes<REAL>> outRec(in[0].size());
	std::vector<DOF_6D_Modes<REAL>> outLig(in[1].size());

		for (int i=0; i < in[0].size(); i++){

				outRec[i]._6D.ang.x = in[0][i]._6D.ang.x;
				outRec[i]._6D.ang.y = in[0][i]._6D.ang.y;
				outRec[i]._6D.ang.z = in[0][i]._6D.ang.z;
				outRec[i]._6D.pos.x = in[0][i]._6D.pos.x;
				outRec[i]._6D.pos.y = in[0][i]._6D.pos.y;
				outRec[i]._6D.pos.z = in[0][i]._6D.pos.z;
		}

		for (int i=0; i < in[1].size(); i++){
				outLig[i]._6D.ang.x = in[1][i]._6D.ang.x;
				outLig[i]._6D.ang.y = in[1][i]._6D.ang.y;
				outLig[i]._6D.ang.z = in[1][i]._6D.ang.z;
				outLig[i]._6D.pos.x = in[1][i]._6D.pos.x;
				outLig[i]._6D.pos.y = in[1][i]._6D.pos.y;
				outLig[i]._6D.pos.z = in[1][i]._6D.pos.z;



				for (int mode=0; mode < Common_Modes::numModesRec; mode++){
					if ( in[1][i].numDofs < Common_Modes::numModesRec || std::isnan(in[1][i].dofs[mode])){
						outLig[i].modesRec[mode] = 0;
					}
					else{
						outLig[i].modesRec[mode] =  in[1][i].dofs[mode];
					}
				}
				for (int mode=0; mode < Common_Modes::numModesLig; mode++){
					if ( in[1][i].numDofs < Common_Modes::numModesLig || std::isnan(in[1][i].dofs[mode])){
						outLig[i].modesLig[mode] = 0;
					}
					else{
						outLig[i].modesLig[mode] =  in[1][i].dofs[mode];
					}
				}


			}
		std::vector<std::vector<DOF_6D_Modes<REAL>>> outVec;
		outVec.push_back(outRec);
		outVec.push_back(outLig);
		return outVec;
}

}// namespace as


#endif /* DOFCONVERTER_H_ */
