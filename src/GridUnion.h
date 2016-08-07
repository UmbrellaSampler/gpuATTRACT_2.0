/*
 * GridUnion.h
 *
 *  Created on: Jul 19, 2016
 *      Author: uwe
 */

#ifndef SRC_GRIDUNION_H_
#define SRC_GRIDUNION_H_

#include <memory>
#include "IntrplGrid.h"
#include "NLGrid.h"
#include "DataItem.h"

namespace as {

template<typename REAL>
class GridUnion : public DataItem {
public:
	/* Constructor */
	GridUnion() {}

	GridUnion(std::shared_ptr<IntrplGrid<REAL>> innerGrid, std::shared_ptr<IntrplGrid<REAL>> outerGrid,
              std::shared_ptr<NLGrid<REAL>> NLGrid) :
            	  inner(innerGrid), outer(outerGrid),
            	  NL(NLGrid) {}

	virtual ~GridUnion() {}

	bool operator ==( const GridUnion& rightSide) const {
			return _tag == rightSide._tag;
	}

	std::string tag() const {
		return _tag;
	}

	void setTag(std::string tag) {
		_tag = tag;
	}

	void freeHost() {
		inner.reset();
		outer.reset();
		NL.reset();
	}

	std::shared_ptr<IntrplGrid<REAL>> inner;
	std::shared_ptr<IntrplGrid<REAL>> outer;
	std::shared_ptr<NLGrid<REAL>> NL;
private:
	std::string _tag;
};

}  // namespace as



#endif /* SRC_GRIDUNION_H_ */
