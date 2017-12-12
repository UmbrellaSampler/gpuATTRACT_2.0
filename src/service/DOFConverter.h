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
	std::vector<std::vector<DOF_6D<REAL>>> outVec(in.size());
	for (int i; i < in.size(); i++){
		for (int k; k < in[i].size(); k++){
			DOF_6D<REAL> out;
			out.ang.x = in[i][k]._6D.ang.x;
			out.ang.y = in[i][k]._6D.ang.y;
			out.ang.z = in[i][k]._6D.ang.z;
			out.pos.x = in[i][k]._6D.pos.x;
			out.pos.y = in[i][k]._6D.pos.y;
			out.pos.z = in[i][k]._6D.pos.z;
			outVec[i].push_back(out);
		}
	}
return outVec;
}


template<typename REAL>
std::vector<std::vector<DOF_6D_Modes<REAL>>> DOFConverter_Modes(std::vector<std::vector<DOF>>  &in) {
	std::vector<std::vector<DOF_6D_Modes<REAL>>> outVec(in.size());

		for (int i; i < in[1].size(); i++){
				DOF_6D_Modes<REAL> outRec;
				outRec._6D.ang.x = in[0][i]._6D.ang.x;
				outRec._6D.ang.y = in[0][i]._6D.ang.y;
				outRec._6D.ang.z = in[0][i]._6D.ang.z;
				outRec._6D.pos.x = in[0][i]._6D.pos.x;
				outRec._6D.pos.y = in[0][i]._6D.pos.y;
				outRec._6D.pos.z = in[0][i]._6D.pos.z;
				outVec[0].push_back(outRec);

				DOF_6D_Modes<REAL> outLig;
				outLig._6D.ang.x = in[1][i]._6D.ang.x;
				outLig._6D.ang.y = in[1][i]._6D.ang.y;
				outLig._6D.ang.z = in[1][i]._6D.ang.z;
				outLig._6D.pos.x = in[1][i]._6D.pos.x;
				outLig._6D.pos.y = in[1][i]._6D.pos.y;
				outLig._6D.pos.z = in[1][i]._6D.pos.z;
				float test=in[1][i]._6D.pos.z;

				for (int mode=0; mode < Common_Modes::numModesRec; mode++){
					if ( in[1][i].numDofs < Common_Modes::numModesRec || std::isnan(in[1][i].dofs[mode])){
						outLig.modesRec[mode] = 0;
					}
					else{
						outLig.modesRec[mode] =  in[1][i].dofs[mode];
					}
				}
				for (int mode=0; mode < Common_Modes::numModesLig; mode++){
					if ( in[1][i].numDofs < Common_Modes::numModesLig || std::isnan(in[1][i].dofs[mode])){
						outLig.modesLig[mode] = 0;
					}
					else{
						outLig.modesLig[mode] =  in[1][i].dofs[mode];
					}
				}

				outVec[1].push_back(outLig);
			}

		return outVec;
}

}// namespace as


#endif /* DOFCONVERTER_H_ */
