#include "AtmosphericProperty1D.h"
#include "units.h"
#include "util.h"
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iostream>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"


NCPA::AtmosphericProperty1D::AtmosphericProperty1D() {
	z_ = NULL;
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units ) {

	z_ = new double[ n_points ];
	z_units_ = altitude_units;
	std::memcpy( z_, altitude_points, n_points*sizeof(double) );
	
	values_ = new double[ n_points ];
	units_ = property_units;
	std::memcpy( values_, property_values, n_points*sizeof(double) );
	n_ = n_points;

	build_splines_();
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( const AtmosphericProperty1D &source ) : VectorWithUnits( source ) {
	this->n_ = source.n_;
	this->z_ = new double[ n_ ];
	std::memcpy( z_, source.z_, n_ * sizeof( double ) );
	this->z_units_ = source.z_units_;

	build_splines_();
}

NCPA::AtmosphericProperty1D::~AtmosphericProperty1D() {
	delete_splines_();
	delete [] z_;
	//delete [] values_;
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_altitude_units() const {
	return z_units_;
}

void NCPA::AtmosphericProperty1D::convert_altitude_units( NCPA::units_t new_units ) {

	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != z_units_) {
		do_units_conversion_( n_, z_, z_units_, new_units );
		build_splines_();
		z_units_ = new_units;
	}
}


void NCPA::AtmosphericProperty1D::convert_units( NCPA::units_t new_units ) {
	//std::cout << "Called AtmosphericProperty1D::convert_units()" << std::endl;
	NCPA::VectorWithUnits::convert_units( new_units );
	build_splines_();
}

void NCPA::AtmosphericProperty1D::build_splines_() {
	delete_splines_();
	// construct spline
	accel_ = gsl_interp_accel_alloc();
#if GSL_MAJOR_VERSION > 1
	spline_ = gsl_spline_alloc( gsl_interp_steffen, n_ );
#else
	spline_ = gsl_spline_alloc( gsl_interp_cspline, n_ );
#endif
	gsl_spline_init( spline_, z_, values_, n_ );
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
	*buffer_units = z_units_;
	std::memcpy( buffer, z_, n_ * sizeof(double) );
}


int NCPA::AtmosphericProperty1D::check_altitude_( double z_req ) const {
	if ( z_req < z_[0] ) {
		return -1;
	} else if ( z_req > z_[ n_ - 1 ] ) {
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
		return values_[ n_ - 1 ];
	} else if ( check < 0 ) {
		return values_[ 0 ];
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
	double z0 = z_[0];
	double z1 = z_[ n_ - 1 ];
	size_t new_nz = (size_t)std::floor( ( z1 - z0 ) / new_dz ) + 1;

	double *new_z = new double[ new_nz ];
	double *new_prop = new double[ new_nz ];
	for (size_t i = 0; i < new_nz; i++) {
		new_z[ i ] = z0 + ((double)i) * new_dz;
		new_prop[ i ] = get( new_z[ i ] );
	}

	// reset everything
	delete [] z_;
	delete [] values_;
	n_ = new_nz;
	z_ = new_z;
	values_ = new_prop;
	build_splines_();
}
