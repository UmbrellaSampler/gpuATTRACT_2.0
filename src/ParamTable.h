#ifndef PARAMTABLE_H_
#define PARAMTABLE_H_

#include <stdexcept>
#include "nativeTypesWrapper.h"

#include "DataItem.h"

namespace as {

/* Potential shape type of the ATTRACT force field */
enum class PotShape {
	_12_6,
	_8_6,
	undefined
};

template<typename REAL>
class ParamTable : public DataItem {
	// Check if REAL is of floating-point type
	using real_t = typename TypeWrapper<REAL>::real_t;
public:

	struct attractFFParams_t {
		real_t rc; // A = rc
		real_t ac; // B = ac
		int ipon;
		real_t rmin2;
		real_t emin;
	};

	using type_t = attractFFParams_t;

	ParamTable() : _paramTable(nullptr), _numTypes(0),
			_shape(PotShape::undefined), _swiOn(0), _swiOff(0) {}

	~ParamTable() {
		delete[] _paramTable;
	}

	unsigned numTypes() const noexcept {
		return _numTypes;
	}

	/* read only access */
	const type_t* table() const noexcept {
		return _paramTable;
	}

	PotShape potShape() const noexcept {
		return _shape;
	}

	real_t swiOn() const noexcept {
		return _swiOn;
	}

	real_t swiOff() const noexcept {
		return _swiOff;
	}

	void setNumTypes(unsigned numTypes) noexcept {
		_numTypes = numTypes;
	}

	void setPotShape(PotShape shape) noexcept {
		_shape = shape;
	}

	void setSwiOn(real_t swiOn) noexcept {
		_swiOn = swiOn;
	}

	void setSwiOff(real_t swiOff) noexcept {
		_swiOff = swiOff;
	}

	const type_t& getParams(const int& typeA, const int& typeB) const noexcept {
		return _paramTable[_numTypes*typeA + typeB];
	}

	type_t* getOrCreateTable() {
		if (_paramTable == nullptr) {
			if (_numTypes == 0) {
				throw std::runtime_error("Error: getOrCreateTypePtr(): the number of types must be set before");
			}
			_paramTable = new ParamTable::type_t[_numTypes * _numTypes];
		}
		return _paramTable;
	}

private:

	type_t* _paramTable;
	unsigned _numTypes; /** number of particle/atom types */

	PotShape _shape; /** potential shape 12:6 or 8:6 or undefined */

	real_t _swiOn;	/** switching potential parameter. Not used at the moment */
	real_t _swiOff;
};

}


#endif /* PARAMTABLE_H_ */
