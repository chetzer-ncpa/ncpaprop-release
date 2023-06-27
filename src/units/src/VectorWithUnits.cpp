#include "units.h"
#include "VectorWithUnits.h"
#include <cstring>
#include <algorithm>


/*
 * Constructors
 */

NCPA::VectorWithUnits::VectorWithUnits() : std::vector<NCPA::ScalarWithUnits>() {}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values,
		units_t property_units )
	: std::vector<NCPA::ScalarWithUnits>(n_points) {
	this->set( n_points, property_values, property_units );
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values,
		const std::string &property_units )
	: VectorWithUnits( n_points, property_values, NCPA::Units::fromString(property_units) ) {
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const NCPA::ScalarWithUnits &singlevalue )
	: std::vector<NCPA::ScalarWithUnits>( n_points, singlevalue ) {}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points,
		const NCPA::ScalarWithUnits *scalarvalues )
	: std::vector<NCPA::ScalarWithUnits>( n_points ) {
	for (size_t i = 0; i < n_points; i++) {
		this->at(i) = scalarvalues[i];
	}
	normalize_units();
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, double singleValue, units_t units )
	: std::vector<NCPA::ScalarWithUnits>( n_points, NCPA::ScalarWithUnits(singleValue, units) ) {}

NCPA::VectorWithUnits::VectorWithUnits( const NCPA::VectorWithUnits &source )
	: std::vector<NCPA::ScalarWithUnits>(source) {
	normalize_units();
}

NCPA::VectorWithUnits::VectorWithUnits( NCPA::VectorWithUnits &&source ) noexcept
	: std::vector<NCPA::ScalarWithUnits>() {
	::swap(*this,source);
}

NCPA::VectorWithUnits::~VectorWithUnits() {
	this->clear();
}

/*
 * Operators
 */

NCPA::VectorWithUnits& NCPA::VectorWithUnits::operator=( NCPA::VectorWithUnits other ) {
	::swap(*this, other);
	return *this;
}


/*
 * Friends
 */

void swap(NCPA::VectorWithUnits &a, NCPA::VectorWithUnits &b) noexcept {
	using std::swap;
	swap(static_cast<std::vector<NCPA::ScalarWithUnits>&>(a),
			static_cast<std::vector<NCPA::ScalarWithUnits>&>(b) );
}


/*
 * Methods
 */
void NCPA::VectorWithUnits::as_array( NCPA::ScalarWithUnits *&buffer, bool normFirst ) {
	if (normFirst) {
		this->normalize_units();
	}
	if (buffer == nullptr) {
		buffer = new NCPA::ScalarWithUnits[ this->size() ];
	}
	size_t i = 0;
	for (auto it = this->cbegin(); it != this->cend(); ++it) {
		buffer[ i++ ] = *it;
	}
}

void NCPA::VectorWithUnits::as_array( double *&buffer, NCPA::units_t &units, bool normFirst ) {
	if (normFirst) {
		this->normalize_units();
	} else if (!this->is_normalized()) {
		throw std::logic_error( "Multiple units present in vector, normalize first!" );
	}
	if (buffer == nullptr) {
		buffer = new double[ this->size() ];
	}
	this->get_values( buffer );
	units = this->get_units();
}


void NCPA::VectorWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw invalid_conversion and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation
	NCPA::units_t oldunits = this->get_units();
	if (new_units != oldunits) {
		NCPA::VectorWithUnits buffer( *this );
		for (auto it = buffer.begin(); it != buffer.end(); ++it) {
			it->convert_units( new_units );
		}
		std::swap(*this,buffer);
	}
}

void NCPA::VectorWithUnits::convert_units( const std::string &new_units ) {
	this->convert_units( NCPA::Units::fromString( new_units ) );
}

void NCPA::VectorWithUnits::fill( const NCPA::ScalarWithUnits &constant ) {
	for (auto it = this->begin(); it != this->end(); ++it) {
		*it = constant;
	}
}

void NCPA::VectorWithUnits::fill( double value, NCPA::units_t units ) {
	this->fill( ScalarWithUnits( value, units ) );
}

void NCPA::VectorWithUnits::fill( double value, const std::string &units ) {
	this->fill( ScalarWithUnits( value, units ) );
}

NCPA::units_t NCPA::VectorWithUnits::get_units( bool normFirst ) {
	if (normFirst) {
		this->normalize_units();
	} else if (!this->is_normalized()) {
		throw std::logic_error( "Multiple units present in vector, normalize first!" );
	}
	return this->begin()->get_units();
}

void NCPA::VectorWithUnits::get_values( size_t &n, double* buffer, bool normFirst ) {
	if (normFirst) {
		this->normalize_units();
	}
	size_t i = 0;
	for (auto it = this->cbegin(); it != this->cend(); ++it) {
		buffer[ i++ ] = it->get();
	}
	n = this->size();
}

void NCPA::VectorWithUnits::get_values( double* buffer, bool normFirst ) {
	size_t n;
	this->get_values( n, buffer, normFirst );
}

bool NCPA::VectorWithUnits::is_normalized() const {
	NCPA::units_t base = this->front().get_units();
	for (auto it = this->cbegin(); it != this->cend(); ++it) {
		NCPA::units_t u = it->get_units();
		if (u != base) {
			return false;
		}
	}
	return true;
}

void NCPA::VectorWithUnits::normalize_units() {
	NCPA::units_t base = this->front().get_units();
	for (auto it = this->begin()+1; it != this->end(); ++it) {
		it->convert_units( base );
	}
}

void NCPA::VectorWithUnits::set( size_t n_points, const double *property_values,
		NCPA::units_t units ) {
	this->resize( n_points );
	for (size_t i = 0; i < n_points; i++) {
		this->at(i) = ScalarWithUnits( property_values[i], units );
	}
}

void NCPA::VectorWithUnits::set( size_t n_points, const double *property_values,
		const std::string &units ) {
	this->set( n_points, property_values, NCPA::Units::fromString( units ) );
}

void NCPA::VectorWithUnits::set( size_t n_points, const ScalarWithUnits *values ) {
	this->resize( n_points );
	for (size_t i = 0; i < n_points; i++) {
		this->at(i) = values[i];
	}
}

void NCPA::VectorWithUnits::set_units( NCPA::units_t new_units ) {
	for (auto it = this->begin(); it != this->end(); ++it) {
		it->set_units( new_units );
	}
}

void NCPA::VectorWithUnits::set_units( const std::string &new_units ) {
	for (auto it = this->begin(); it != this->end(); ++it) {
		it->set_units( new_units );
	}
}



