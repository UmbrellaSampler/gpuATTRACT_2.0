#ifndef PARAMTABLE_H_
#define PARAMTABLE_H_

#include <stdexcept>

namespace as {

template<typename REAL>
class AttrParamTable {
public:

	template<typename _REAL>
	struct attractFFParams_t {
		_REAL rc; // A = rc
		_REAL ac; // B = ac
		int ipon;
		_REAL rmin2;
		_REAL emin;
	};

	/* Potential shape type of the ATTRACT force field */
	enum PotShape {
		_12_6,
		_8_6,
		undefined
	};

	using type = attractFFParams_t<REAL>;

	/* Constructor */
	AttrParamTable() : _paramTable(nullptr), _numTypes(0),
			_shape(undefined), _swiOn(0), _swiOff(0) {}



	/* Destructor */
	~AttrParamTable() {
		delete[] _paramTable;
	}

	/***************
	* G E T T E R
	***************/
	unsigned numTypes() const noexcept {
		return _numTypes;
	}

	/* read only access */
	const type* table() const noexcept {
		return _paramTable;
	}

	PotShape potShape() const noexcept {
		return _shape;
	}

	REAL swiOn() const noexcept {
		return _swiOn;
	}

	REAL swiOff() const noexcept {
		return _swiOff;
	}

	void setNumTypes(unsigned numTypes) noexcept {
		_numTypes = numTypes;
	}

	void setPotShape(PotShape shape) noexcept {
		_shape = shape;
	}

	void setSwiOn(REAL swiOn) noexcept {
		_swiOn = swiOn;
	}

	void setSwiOff(REAL swiOff) noexcept {
		_swiOff = swiOff;
	}

	/****************************
	 * public member functions
	 ****************************/
	inline const type& getParams(const int& typeA, const int& typeB) const noexcept {
		return _paramTable[_numTypes*typeA + typeB];
	}

	/*
	 * Read and write access.
	 * Should be used for initialization
	 */
	type* getOrCreateTable() {
		if (_paramTable == nullptr) {
			if (_numTypes == 0) {
				throw std::runtime_error("Error: getOrCreateTypePtr(): the number of types must be set before");
			}
			_paramTable = new AttrParamTable::type[_numTypes * _numTypes];
		}
		return _paramTable;
	}

private:

	type* _paramTable;
	unsigned _numTypes; /** number of particle/atom types */

	PotShape _shape; /** potential shape 12:6 or 8:6 or undefined */

	REAL _swiOn;	/** switching potential parameter. Not used at the moment */
	REAL _swiOff;
};

}


#endif /* PARAMTABLE_H_ */
