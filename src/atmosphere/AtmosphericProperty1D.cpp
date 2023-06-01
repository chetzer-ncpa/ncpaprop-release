#include "AtmosphericProperty1D.h"
#include "units.h"
#include "util.h"
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"


NCPA::AtmosphericProperty1D::AtmosphericProperty1D() {
//	z_ = NULL;
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units ) {

//	z_ = new double[ n_points ];
//	z_units_ = altitude_units;
//	std::memcpy( z_, altitude_points, n_points*sizeof(double) );
//	z_.( n_points, altitude_points, altitude_units );
//
//	values_ = new double[ n_points ];
//	units_ = property_units;
//	std::memcpy( values_, property_values, n_points*sizeof(double) );
//	n_ = n_points;
	this->first.set( n_points, altitude_points, altitude_units );
	this->second.set( n_points, property_values, property_units );
	build_splines_();
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( const AtmosphericProperty1D &source ) {
	this->first = source.first;
	this->second = source.second;
//	this->n_ = source.n_;
//	this->z_ = new double[ n_ ];
//	std::memcpy( z_, source.z_, n_ * sizeof( double ) );
//	this->z_units_ = source.z_units_;

	build_splines_();
}

NCPA::AtmosphericProperty1D::~AtmosphericProperty1D() {
	delete_splines_();
//	delete [] z_;
	//delete [] values_;
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_altitude_units() const {
//	return z_units_;
	return this->first.get_units();
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_units() const {
//	return z_units_;
	return this->second.get_units();
}

void NCPA::AtmosphericProperty1D::convert_altitude_units( NCPA::units_t new_units ) {

	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
//	if (new_units != z_units_) {
//		do_units_conversion_( n_, z_, z_units_, new_units );
//		build_splines_();
//		z_units_ = new_units;
//	}
	if (new_units != this->first.get_units()) {
		this->first.convert_units( new_units );
		build_splines_();
	}
}


void NCPA::AtmosphericProperty1D::convert_units( NCPA::units_t new_units ) {
	//std::cout << "Called AtmosphericProperty1D::convert_units()" << std::endl;
//	NCPA::VectorWithUnits::convert_units( new_units );
	this->second.convert_units( new_units );
	build_splines_();
}

void NCPA::AtmosphericProperty1D::build_splines_() {
	delete_splines_();
	// construct spline
	size_t n = this->first.size();
	accel_ = gsl_interp_accel_alloc();
	spline_ = gsl_spline_alloc( ATMOSPHERIC_INTERPOLATION_TYPE, n );
	double *z_buffer = new double[ n ];
	double *x_buffer = new double[ n ];
	this->first.get_vector( z_buffer );
	this->second.get_vector( x_buffer );
	gsl_spline_init( spline_, z_buffer, x_buffer, n );
	delete [] z_buffer;
	delete [] x_buffer;
}

void NCPA::AtmosphericProperty1D::delete_splines_() {
	if (spline_ != NULL) {
		gsl_spline_free( spline_ );
		spline_ = NULL;
	}
	if (accel_ != NULL) {
		gsl_interp_accel_free( accel_ );
		accel_ = NULL;
	}
}

void NCPA::AtmosphericProperty1D::get_altitude_vector( double *buffer, units_t *buffer_units ) const {
//	*buffer_units = z_units_;
//	std::memcpy( buffer, z_, n_ * sizeof(double) );
	this->first.get_vector( buffer, buffer_units );
}

void NCPA::AtmosphericProperty1D::get_altitude_vector( double *buffer ) const {
//	*buffer_units = z_units_;
//	std::memcpy( buffer, z_, n_ * sizeof(double) );
	this->first.get_vector( buffer );
}


void NCPA::AtmosphericProperty1D::get_vector( double *buffer, units_t *buffer_units ) const {
//	*buffer_units = z_units_;
//	std::memcpy( buffer, z_, n_ * sizeof(double) );
	this->second.get_vector( buffer, buffer_units );
}

void NCPA::AtmosphericProperty1D::get_vector( double *buffer ) const {
//	*buffer_units = z_units_;
//	std::memcpy( buffer, z_, n_ * sizeof(double) );
	this->second.get_vector( buffer );
}

int NCPA::AtmosphericProperty1D::check_altitude_( double z_req ) const {
	if ( z_req < this->first.front() ) {
		return -1;
	} else if ( z_req > this->first.back() ) {
		return 1;
	} else {
		return 0;
	}
	// if ( z_req < z_[0] || z_req > z_[ n_ - 1 ] ) {
	// 	std::ostringstream oss;
	// 	oss << "Requested altitude " << z_req << " " << NCPA::Units::toStr( z_units_ ) << " outside profile bounds.";
	// 	throw std::range_error( oss.str() );
	// }
}

double NCPA::AtmosphericProperty1D::get( double z_req ) const {
	// check_altitude_( z_req );
	int check = check_altitude_( z_req );
	if ( check > 0 ) {
		return this->second.back();
//		return values_[ n_ - 1 ];
	} else if ( check < 0 ) {
		return this->second.front();
//		return values_[ 0 ];
	} else {
		return gsl_spline_eval( spline_, z_req, accel_ );
	}
}

double NCPA::AtmosphericProperty1D::get_first_derivative( double z_req ) const {
	// check_altitude_( z_req );
	int check = check_altitude_( z_req );
	if (check == 0) {
		return gsl_spline_eval_deriv( spline_, z_req, accel_ );
	} else {
		return 0;
	}
}

double NCPA::AtmosphericProperty1D::get_second_derivative( double z_req ) const {
	// check_altitude_( z_req );
	int check = check_altitude_( z_req );
	if (check == 0) {
		return gsl_spline_eval_deriv2( spline_, z_req, accel_ );
	} else {
		return 0;
	}
}

void NCPA::AtmosphericProperty1D::resample( double new_dz ) {

	// get new altitude vector
//	double z0 = z_[0];
//	double z1 = z_[ n_ - 1 ];
	double z0 = this->first.front();
	double z1 = this->first.back();
	size_t new_nz = (size_t)std::floor( ( z1 - z0 ) / new_dz ) + 1;

	double *new_z = new double[ new_nz ];
	double *new_prop = new double[ new_nz ];
	for (size_t i = 0; i < new_nz; i++) {
		new_z[ i ] = z0 + ((double)i) * new_dz;
		new_prop[ i ] = this->get( new_z[ i ] );
	}

	// reset everything
//	delete [] z_;
//	delete [] values_;
//	n_ = new_nz;
//	z_ = new_z;
//	values_ = new_prop;
	this->first.set_values( new_nz, new_z );
	this->second.set_values( new_nz, new_prop );
	build_splines_();
}
