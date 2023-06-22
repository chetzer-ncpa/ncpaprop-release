#include "units.h"
#include "VectorWithUnits.h"
#include <cstring>
#include <algorithm>



NCPA::VectorWithUnits::VectorWithUnits() : std::vector<double>(), units_{ NCPA::units_t::NONE } {}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values, units_t property_units )
	: std::vector<double>(n_points) {
	this->set( n_points, property_values, property_units );
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values, const std::string &property_units )
	: VectorWithUnits( n_points, property_values, NCPA::Units::fromString(property_units) ) {
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const NCPA::ScalarWithUnits &singlevalue )
	: std::vector<double>( n_points, singlevalue.get() ), units_{singlevalue.get_units()} {}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const NCPA::ScalarWithUnits *scalarvalues )
	: std::vector<double>( n_points ) {
	set_units( scalarvalues[0].get_units() );
	for (size_t i = 0; i < this->size(); i++) {
		this->at(i) = NCPA::Units::convert( scalarvalues[i].get(), scalarvalues[i].get_units(), units_ );
	}
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, double singleValue, units_t units )
	: std::vector<double>( n_points, singleValue ), units_{units} {}

NCPA::VectorWithUnits::VectorWithUnits( const NCPA::VectorWithUnits &source ) : std::vector<double>(source) {
	this->set_units( source.get_units() );
}

NCPA::VectorWithUnits::VectorWithUnits( NCPA::VectorWithUnits &&source ) noexcept : std::vector<double>() {
	::swap(*this,source);
}

NCPA::VectorWithUnits::~VectorWithUnits() {
	this->clear();
}

void swap(NCPA::VectorWithUnits &a, NCPA::VectorWithUnits &b) noexcept {
	using std::swap;
	swap(static_cast<std::vector<double>&>(a), static_cast<std::vector<double>&>(b) );
	swap(a.units_,b.units_);
}

void NCPA::VectorWithUnits::set( size_t n_points, const double *property_values, NCPA::units_t units ) {
	this->set_units( units );
	this->set_values( n_points, property_values );
}

void NCPA::VectorWithUnits::set_values( size_t n_points, const double *property_values ) {
//	this->clear();
//	this->resize( n_points );
//	std::copy(property_values, property_values + n_points, this->begin() );
	this->assign(property_values, property_values + n_points);
}

void NCPA::VectorWithUnits::fill( size_t n_points, double value ) {
	this->resize( n_points );
	this->fill( value );
}

void NCPA::VectorWithUnits::fill( double value ) {
	for (auto it = this->begin(); it != this->end(); ++it) {
		*it = value;
	}
}

NCPA::units_t NCPA::VectorWithUnits::get_units() const {
	return units_;
}

void NCPA::VectorWithUnits::set_units( NCPA::units_t new_units ) {
	units_ = new_units;
}

void NCPA::VectorWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw invalid_conversion and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation
	if (new_units != units_) {
		std::vector<double> buffer( this->begin(), this->end() );
		for (std::vector<double>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
			*it = NCPA::Units::convert( *it, units_, new_units );
		}
		this->assign(buffer.begin(),buffer.end());
		this->set_units( new_units );
	}
}

void NCPA::VectorWithUnits::convert_units( const std::string &new_units ) {
	this->convert_units( NCPA::Units::fromString( new_units ) );
}


void NCPA::VectorWithUnits::get_vector( double *buffer, units_t &buffer_units ) const {
	buffer_units = units_;
	std::copy(this->begin(), this->end(), buffer);
}

void NCPA::VectorWithUnits::get_vector( double *buffer ) const {
	std::copy(this->begin(), this->end(), buffer);
}


NCPA::VectorWithUnits& NCPA::VectorWithUnits::operator=( NCPA::VectorWithUnits other ) {
	::swap(*this, other);
	return *this;
}
