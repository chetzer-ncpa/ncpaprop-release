#include "units.h"
#include "VectorWithUnits.h"
#include <cstring>
#include <algorithm>



NCPA::VectorWithUnits::VectorWithUnits() : units_{ NCPA::UNITS_NONE } {}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values, units_t property_units ) {
	this->set( n_points, property_values, property_units );
}



NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values, const std::string &property_units )
	: VectorWithUnits( n_points, property_values, NCPA::Units::fromString(property_units) ) {
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const NCPA::ScalarWithUnits &singlevalue ) {
	this->set_units( singlevalue.get_units() );
	this->fill( n_points, singlevalue.get() );
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const NCPA::ScalarWithUnits *scalarvalues ) {
	this->clear();
	this->resize( n_points );
	set_units( scalarvalues[0].get_units() );
	for (size_t i = 0; i < this->size(); i++) {
		this->at(i) = NCPA::Units::convert( scalarvalues[i].get(), scalarvalues[i].get_units(), units_ );
	}
}

NCPA::VectorWithUnits::VectorWithUnits( const NCPA::VectorWithUnits &source ) {
	this->set_units( source.get_units() );
	this->assign(source.begin(),source.end());
}

NCPA::VectorWithUnits::~VectorWithUnits() {
	this->clear();
}

void NCPA::VectorWithUnits::set( size_t n_points, const double *property_values, NCPA::units_t units ) {
	this->set_units( units );
	this->set_values( n_points, property_values );
}

void NCPA::VectorWithUnits::set_values( size_t n_points, const double *property_values ) {
	this->clear();
	this->resize( n_points );
	std::copy(property_values, property_values + n_points, this->begin() );
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
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation
	if (new_units != units_) {
		std::vector<double> buffer( this->begin(), this->end() );
		for (std::vector<double>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
			*it = NCPA::Units::convert( *it, units_, new_units );
		}
		this->assign(buffer.begin(),buffer.end());
		this->set_units( new_units );
	}
//	if (new_units != units_) {
//		do_units_conversion_( n_, values_, units_, new_units );
//		units_ = new_units;
//	}
}

void NCPA::VectorWithUnits::convert_units( const std::string &new_units ) {
	this->convert_units( NCPA::Units::fromString( new_units ) );
}

//void NCPA::VectorWithUnits::do_units_conversion_( size_t n_points, double *inplace,
//			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {
//
//	// try to convert
//	double *units_buffer = new double[ n_points ];
//	std::memset( units_buffer, 0, n_points * sizeof( double ) );
//
//	// throws out_of_range if conversion is undefined
//	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );
//
//	// successful, so record the units change
//	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
//	delete [] units_buffer;
//}


//size_t NCPA::VectorWithUnits::size() const {
//	return n_;
//}

void NCPA::VectorWithUnits::get_vector( double *buffer, units_t *buffer_units ) const {
	*buffer_units = units_;
	std::copy(this->begin(), this->end(), buffer);
//	std::memcpy( buffer, values_, n_ * sizeof( double ) );
}

void NCPA::VectorWithUnits::get_vector( double *buffer ) const {
	std::copy(this->begin(), this->end(), buffer);
//	std::memcpy( buffer, values_, n_ * sizeof( double ) );
}

//double NCPA::VectorWithUnits::operator[]( size_t i ) const {
//	return values_[ i ];
//}
//
//double &NCPA::VectorWithUnits::operator[]( size_t i ) {
//	return values_[ i ];
//}

NCPA::VectorWithUnits& NCPA::VectorWithUnits::operator=( const NCPA::VectorWithUnits& other ) {
	if (this == &other) return *this;

	double *buffer = new double[ other.size() ];
	NCPA::units_t otherunits;
	other.get_vector( buffer, &otherunits );
	this->set( other.size(), buffer, otherunits );
	delete [] buffer;
	return *this;
}
