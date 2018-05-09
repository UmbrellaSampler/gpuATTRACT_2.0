#ifndef RETRANSFORMDOFS_TPP_
#define RETRANSFORMDOFS_TPP_


#include <vector>
#include <Types_6D.h>
#include <Types_6D_Modes.h>

namespace as{


template <typename REAL>
void transformDOF(	bool input_centered_receptor, bool input_centered_ligand,
					bool output_centered_receptor, bool output_centered_ligand,
					Vec3<double> pivot_receptor,Vec3<double> pivot_ligand,
					std::vector<DOF_6D<REAL>> &dofs){
	for( auto & dof: dofs){
		if (!input_centered_ligand) {
			dof.pos.x -= pivot_ligand.x;
			dof.pos.y -= pivot_ligand.y;
			dof.pos.z -= pivot_ligand.z;
			}
		else{
			dof.pos.x -= pivot_ligand.x;
			dof.pos.y -= pivot_ligand.y;
			dof.pos.z -= pivot_ligand.z;
		}
		/* shift ligand dofs by receptor pivot */
		if (!input_centered_receptor) {
			dof.pos.x += pivot_receptor.x;
			dof.pos.y += pivot_receptor.y;
			dof.pos.z += pivot_receptor.z;
		}

		if (output_centered_ligand) {
			dof.pos.x -= pivot_ligand.x;
			dof.pos.y -= pivot_ligand.y;
			dof.pos.z -= pivot_ligand.z;
		}

		/* shift ligand dofs by receptor pivot */
		if (output_centered_receptor) {
			dof.pos.x += pivot_receptor.x;
			dof.pos.y += pivot_receptor.y;
			dof.pos.z += pivot_receptor.z;
		}
	}
}

template <typename REAL>
void transformDOF(	bool input_centered_receptor, bool input_centered_ligand,
					bool output_centered_receptor, bool output_centered_ligand,
					Vec3<double> pivot_receptor,Vec3<double> pivot_ligand,
					std::vector<DOF_6D_Modes<REAL>> &dofs){
	for( auto & dof: dofs){
		if (!input_centered_ligand) {
			dof._6D.pos.x -= pivot_ligand.x;
			dof._6D.pos.y -= pivot_ligand.y;
			dof._6D.pos.z -= pivot_ligand.z;
			}
		else {
			dof._6D.pos.x -= pivot_ligand.x;
			dof._6D.pos.y -= pivot_ligand.y;
			dof._6D.pos.z -= pivot_ligand.z;
		}

		/* shift ligand dofs by receptor pivot */
		if (!input_centered_receptor) {
			dof._6D.pos.x += pivot_receptor.x;
			dof._6D.pos.y += pivot_receptor.y;
			dof._6D.pos.z += pivot_receptor.z;
		}

		if (output_centered_ligand) {
			dof._6D.pos.x -= pivot_ligand.x;
			dof._6D.pos.y -= pivot_ligand.y;
			dof._6D.pos.z -= pivot_ligand.z;
		}

		/* shift ligand dofs by receptor pivot */
		if (output_centered_receptor) {
			dof._6D.pos.x += pivot_receptor.x;
			dof._6D.pos.y += pivot_receptor.y;
			dof._6D.pos.z += pivot_receptor.z;
		}
	}
}

}
#endif
