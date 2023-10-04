#include "Interpolator1D.h"


// this will only compile if gsl/gsl_spline.h is available
#ifdef HAVE_GSL_INTERPOLATION_LIBRARY

#include "gsl/gsl_spline.h"
#include "Interpolator1D.h"
#include "GSLInterpolator1D.h"
#include "GSLPeriodicInterpolator1D.h"
#include "NCPACommon.h"
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>

NCPA::GSLPeriodicInterpolator1D::GSLPeriodicInterpolator1D(
		const gsl_interp_type *interp_type, NCPA::interpolator1d_t ncpa_type )
		: NCPA::GSLInterpolator1D(interp_type, ncpa_type) {}

NCPA::GSLPeriodicInterpolator1D::GSLPeriodicInterpolator1D(
		const GSLPeriodicInterpolator1D &other )
		: NCPA::GSLInterpolator1D( other ) {}

NCPA::GSLPeriodicInterpolator1D::GSLPeriodicInterpolator1D(
		GSLPeriodicInterpolator1D &&other )
	: NCPA::GSLInterpolator1D() {
	::swap(*this,other);
}

NCPA::GSLPeriodicInterpolator1D::~GSLPeriodicInterpolator1D() {
	this->free();
}

void swap( NCPA::GSLPeriodicInterpolator1D &a, NCPA::GSLPeriodicInterpolator1D &b ) {
	using std::swap;
	::swap( static_cast<NCPA::GSLInterpolator1D&>(a),
			static_cast<NCPA::GSLInterpolator1D&>(b) );
}

double NCPA::GSLPeriodicInterpolator1D::f( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	double xrange = max_x - min_x;

	if (is_extrapolating()) {
		while (x < min_x) {
			x += xrange;
		}
		while (x > max_x) {
			x -= xrange;
		}
		return interpolate_f_( x );
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
}

double NCPA::GSLPeriodicInterpolator1D::df( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	double xrange = max_x - min_x;

	if (is_extrapolating()) {
		while (x < min_x) {
			x += xrange;
		}
		while (x > max_x) {
			x -= xrange;
		}
		return interpolate_df_( x );
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
}

double NCPA::GSLPeriodicInterpolator1D::df( size_t n, double x ) {
	if (n == 0) {
		return this->f( x );
	} else if (n == 1) {
		return this->df( x );
	}

	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	double xrange = max_x - min_x;

	if (is_extrapolating()) {
		while (x < min_x) {
			x += xrange;
		}
		while (x > max_x) {
			x -= xrange;
		}
//		return interpolate_df_( n, x );
		if (n == 2) {
			return interpolate_d2f_( x );
		} else if (n == 3) {
			return interpolate_d3f_( x );
		} else {
			throw std::out_of_range( "Derivative level must be in [0,3]" );
		}
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
}

std::complex<double> NCPA::GSLPeriodicInterpolator1D::cf( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	double xrange = max_x - min_x;

	if (is_extrapolating()) {
		while (x < min_x) {
			x += xrange;
		}
		while (x > max_x) {
			x -= xrange;
		}
		return interpolate_cf_( x );
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
}

std::complex<double> NCPA::GSLPeriodicInterpolator1D::cdf( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	double xrange = max_x - min_x;

	if (is_extrapolating()) {
		while (x < min_x) {
			x += xrange;
		}
		while (x > max_x) {
			x -= xrange;
		}
		return interpolate_cdf_( x );
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
}

std::complex<double> NCPA::GSLPeriodicInterpolator1D::cdf( size_t n, double x ) {
	if (n == 0) {
		return this->f( x );
	} else if (n == 1) {
		return this->df( x );
	}

	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	double xrange = max_x - min_x;

	if (is_extrapolating()) {
		while (x < min_x) {
			x += xrange;
		}
		while (x > max_x) {
			x -= xrange;
		}
		if (n == 2) {
			return interpolate_cd2f_( x );
		} else if (n == 3) {
			return interpolate_cd3f_( x );
		} else {
			throw std::out_of_range( "Derivative level must be in [0,3]" );
		}

	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
}



#endif
