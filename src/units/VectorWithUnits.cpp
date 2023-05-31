#include "units.h"
#include "VectorWithUnits.h"
#include <cstring>



NCPA::VectorWithUnits::VectorWithUnits() : n_{ 0 }, values_{ NULL }, units_{ NCPA::UNITS_NONE } {}


NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values, units_t property_units )
	: n_{ n_points }, units_{ property_units } {
	values_ = new double[ n_ ];
	std::memcpy( values_, property_values, n_points*sizeof(double) );
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const double *property_values, const std::string &property_units )
	: n_{ n_points }, units_{ NCPA::Units::fromString(property_units) } {
	values_ = new double[ n_ ];
	std::memcpy( values_, property_values, n_points*sizeof(double) );
}

NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, const NCPA::ScalarWithUnits *scalarvalues ) {
	n_ = n_points;
	units_ = scalarvalues[0].get_units();
	values_ = new double[ n_ ];
	for (size_t i = 0; i < n_; i++) {
		values_[i] = NCPA::Units::convert( scalar_values[i].get(), scalar_values[i].get_units, units_ );
	}
}

NCPA::VectorWithUnits::VectorWithUnits( const NCPA::VectorWithUnits &source )
	: n_{ source.n_ }, units_{ source.units_ } {
	//n_ = source.n_;
	//units_ = source.units_;
	values_ = new double[ n_ ];
	std::memcpy( values_, source.values_, n_ * sizeof( double ) );
}

NCPA::VectorWithUnits::~VectorWithUnits() {
	if (values_ != NULL) {
		delete [] values_;
	}
}

NCPA::units_t NCPA::VectorWithUnits::get_units() const {
	return units_;
}

void NCPA::VectorWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation
	if (new_units != units_) {
		do_units_conversion_( n_, values_, units_, new_units );
		units_ = new_units;
	}
}

void NCPA::VectorWithUnits::convert_units( const std::string &new_units ) {
	this->convert_units( NCPA::Units::fromString( new_units ) );
}

void NCPA::VectorWithUnits::do_units_conversion_( size_t n_points, double *inplace,
			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double *units_buffer = new double[ n_points ];
	std::memset( units_buffer, 0, n_points * sizeof( double ) );

	// throws out_of_range if conversion is undefined
	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );

	// successful, so record the units change
	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
	delete [] units_buffer;
}


size_t NCPA::VectorWithUnits::size() const {
	return n_;
}

void NCPA::VectorWithUnits::get_vector( double *buffer, units_t *buffer_units ) const {
	*buffer_units = units_;
	std::memcpy( buffer, values_, n_ * sizeof( double ) );
}

void NCPA::VectorWithUnits::get_vector( double *buffer ) const {
	std::memcpy( buffer, values_, n_ * sizeof( double ) );
}

double NCPA::VectorWithUnits::operator[]( size_t i ) const {
	return values_[ i ];
}

double &NCPA::VectorWithUnits::operator[]( size_t i ) {
	return values_[ i ];
}
