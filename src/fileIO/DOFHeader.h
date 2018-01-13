/*
 * DOFHeader.h
 *
 *  Created on: Dec 7, 2017
 *      Author: uwe
 */

#ifndef SRC_FILEIO_DOFHEADER_H_
#define SRC_FILEIO_DOFHEADER_H_

namespace as {

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

} // namespace as



#endif /* SRC_FILEIO_DOFHEADER_H_ */
