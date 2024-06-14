#include "AtmosphericProperty1D.h"
#include "units.h"
#include "util.h"
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iostream>
#include <cassert>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"

NCPA::AtmosphericProperty1D::AtmosphericProperty1D() {
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points,
		double *altitude_points,
		units_t altitude_units,
		double *property_values,
		units_t property_units ) :
		NCPA::AtmosphericProperty1D::AtmosphericProperty1D() {

	z_ = new double[ n_points ];
	z_units_ = altitude_units;
	std::memcpy( z_, altitude_points, n_points * sizeof(double) );

	values_ = new double[ n_points ];
	units_ = property_units;
	std::memcpy( values_, property_values, n_points * sizeof(double) );
	n_ = n_points;

	build_splines_();
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( const AtmosphericProperty1D &source ) :
		VectorWithUnits( source ) {
	this->n_ = source.n_;
	this->z_ = new double[ n_ ];
	std::memcpy( z_, source.z_, n_ * sizeof(double) );
	this->z_units_ = source.z_units_;

	build_splines_();
}

NCPA::AtmosphericProperty1D::~AtmosphericProperty1D() {
	delete_splines_();
	if (z_ != nullptr) {
		delete[] z_;
	}
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
	if (spline_ != nullptr) {
		gsl_spline_free( spline_ );
		spline_ = nullptr;
	}
	if (accel_ != nullptr) {
		gsl_interp_accel_free( accel_ );
		accel_ = nullptr;
	}
}

void NCPA::AtmosphericProperty1D::get_altitude_vector( double *buffer,
		units_t *buffer_units ) const {
	*buffer_units = z_units_;
	std::memcpy( buffer, z_, n_ * sizeof(double) );
}

int NCPA::AtmosphericProperty1D::check_altitude_( double z_req ) const {
	if (z_req < z_[ 0 ]) {
		return -1;
	} else if (z_req > z_[ n_ - 1 ]) {
		return 1;
	} else {
		return 0;
	}
}

double NCPA::AtmosphericProperty1D::get( double z_req ) const {
	int check = check_altitude_( z_req );
	if (check > 0) {
		return values_[ n_ - 1 ];
	} else if (check < 0) {
		return values_[ 0 ];
	} else {
		return gsl_spline_eval( spline_, z_req, accel_ );
	}
}

double NCPA::AtmosphericProperty1D::get_first_derivative( double z_req ) const {
	int check = check_altitude_( z_req );
	if (check > 0) {
		return 0.0;
	} else if (check < 0) {
		return 0.0;
	} else {
		return gsl_spline_eval_deriv( spline_, z_req, accel_ );
	}
}

double NCPA::AtmosphericProperty1D::get_second_derivative( double z_req ) const {
	int check = check_altitude_( z_req );
	if (check > 0) {
		return 0.0;
	} else if (check < 0) {
		return 0.0;
	} else {
		return gsl_spline_eval_deriv2( spline_, z_req, accel_ );
	}
}

void NCPA::AtmosphericProperty1D::resample( double new_dz ) {

	// get new altitude vector
	double z0 = z_[ 0 ];
	double z1 = z_[ n_ - 1 ];
	size_t new_nz = (size_t)std::floor( (z1 - z0) / new_dz ) + 1;

	double *new_z = new double[ new_nz ];
	double *new_prop = new double[ new_nz ];
	for (size_t i = 0; i < new_nz; i++) {
		new_z[ i ] = z0 + ((double)i) * new_dz;
		new_prop[ i ] = get( new_z[ i ] );
	}

	// reset everything
	delete[] z_;
	delete[] values_;
	n_ = new_nz;
	z_ = new_z;
	values_ = new_prop;
	build_splines_();
}

void NCPA::AtmosphericProperty1D::add_point( double new_z, double new_f, bool replace, double tolerance ) {
	this->add_points( 1, &new_z, &new_f, replace, tolerance );
}


void NCPA::AtmosphericProperty1D::add_points( size_t n_new_points,
		const double *z_points,
		const double *f_points,
		bool replace,
		double tolerance ) {
	std::vector<double> new_z, new_f;
	new_z.reserve( n_ + n_new_points );
	new_f.reserve( n_ + n_new_points );

	size_t old_i = 0, new_i = 0;
	while (old_i < n_ || new_i < n_new_points) {
		// case 1: used up all the old points
		if (old_i == n_) {
			new_z.push_back( z_points[ new_i ] );
			new_f.push_back( f_points[ new_i ] );
			new_i++;

			// case 2: used up all the new points
		} else if (new_i == n_new_points) {
			new_z.push_back( z_[ old_i ] );
			new_f.push_back( values_[ old_i ] );
			old_i++;

			// case 3: compare the next points in each vector
		} else {
			double diff = std::fabs( z_[ old_i ] - z_points[ new_i ] );

			// case 3a: repeated altitude value
			if (diff <= tolerance) {
				if (replace) {
					new_z.push_back( z_points[ new_i ] );
					new_f.push_back( f_points[ new_i ] );
					new_i++;
					old_i++;  // cause we're replacing it
				} else {
					std::ostringstream oss;
					oss << "New z value of " << z_points[ new_i ]
						<< " already exists in property vector (" << z_[ old_i ]
						<< " within " << tolerance << ")";
					throw std::range_error( oss.str() );
				}
			} else {
				// case 3b: next z value is in the old altitude vector
				if (z_[ old_i ] < z_points[ new_i ]) {
					new_z.push_back( z_[ old_i ] );
					new_f.push_back( values_[ old_i ] );
					old_i++;

					// case 3c: next z value is in the new altitude vector
				} else {
					new_z.push_back( z_points[ new_i ] );
					new_f.push_back( f_points[ new_i ] );
					new_i++;
				}
			}
		}
	} // end while loop
	assert( new_z.size() == new_f.size() );

	// now reset the internal vectors
	delete[] z_;
	delete[] values_;
	z_ = NCPA::zeros<double>( new_z.size() );
	values_ = NCPA::zeros<double>( new_z.size() );
	std::copy( new_z.begin(), new_z.end(), z_ );
	std::copy( new_f.begin(), new_f.end(), values_ );
	build_splines_();
}

