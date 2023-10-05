#include "Interpolator1D.h"
#include "NCPACommon.h"

#include "../include/NCPALinearInterpolator1D.h"
#include "../include/NCPANearestNeighborInterpolator1D.h"

#ifdef HAVE_LANL_INTERPOLATION_LIBRARY
#include "LANLInterpolation.h"
#include "LANLLinearInterpolator1D.h"
#include "LANLNaturalCubicSplineInterpolator1D.h"
#endif

#ifdef HAVE_GSL_INTERPOLATION_LIBRARY
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"
#include "GSLInterpolator1D.h"
#include "GSLPeriodicInterpolator1D.h"
#endif

#include <complex>
#include <stdexcept>
#include <sstream>
#include <iostream>

std::string NCPA::Interpolator1D::as_string( NCPA::interpolator1d_t t ) {
	switch (t) {
	case NCPA::interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR:
		return "NCPA Nearest-neighbor selector";
		break;
	case NCPA::interpolator1d_t::NCPA_1D_LINEAR:
		return "NCPA linear interpolator";
		break;
	case NCPA::interpolator1d_t::LANL_1D_LINEAR:
		return "LANL linear interpolator";
		break;
	case NCPA::interpolator1d_t::LANL_1D_NATURAL_CUBIC:
		return "LANL natural cubic spline interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_LINEAR:
		return "GNU Scientific Library linear interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_POLYNOMIAL:
		return "GNU Scientific Library polynomial interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_CUBIC:
		return "GNU Scientific Library cubic spline interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_CUBIC_PERIODIC:
		return "GNU Scientific Library periodic cubic spline interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA:
		return "GNU Scientific Library Akima spline interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA_PERIODIC:
		return "GNU Scientific Library periodic Akima spline interpolator";
		break;
	case NCPA::interpolator1d_t::GSL_1D_STEFFEN:
		return "GNU Scientific Library Steffen spline interpolator";
		break;
	default:
		throw std::out_of_range("Unrecognized or unimplemented interpolator requested.");
	}
}

bool NCPA::Interpolator1D::can_build( NCPA::interpolator1d_t t ) {
	switch (t) {
	case NCPA::interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR:
	case NCPA::interpolator1d_t::NCPA_1D_LINEAR:
		return true;
		break;
	case NCPA::interpolator1d_t::LANL_1D_LINEAR:
	case NCPA::interpolator1d_t::LANL_1D_NATURAL_CUBIC:
#ifdef HAVE_LANL_INTERPOLATION_LIBRARY
		return true;
#else
		return false;
#endif
		break;
	case NCPA::interpolator1d_t::GSL_1D_LINEAR:
	case NCPA::interpolator1d_t::GSL_1D_POLYNOMIAL:
	case NCPA::interpolator1d_t::GSL_1D_CUBIC:
	case NCPA::interpolator1d_t::GSL_1D_CUBIC_PERIODIC:
	case NCPA::interpolator1d_t::GSL_1D_AKIMA:
	case NCPA::interpolator1d_t::GSL_1D_AKIMA_PERIODIC:
#ifdef HAVE_GSL_INTERPOLATION_LIBRARY
		return true;
#else
		return false;
#endif
		break;

	case NCPA::interpolator1d_t::GSL_1D_STEFFEN:
#ifdef HAVE_GSL_INTERPOLATION_LIBRARY

#if GSL_MAJOR_VERSION > 1
		return true;
#else
		return false;
#endif

#else
		return false;
#endif
	default:
		return false;
	}
}

NCPA::Interpolator1D *NCPA::Interpolator1D::build( NCPA::interpolator1d_t t ) {
	switch (t) {
	case NCPA::interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR:
		return static_cast<NCPA::Interpolator1D*>(new NCPA::NCPANearestNeighborInterpolator1D());
		break;
	case NCPA::interpolator1d_t::NCPA_1D_LINEAR:
		return static_cast<NCPA::Interpolator1D*>(new NCPA::NCPALinearInterpolator1D());
		break;
#ifdef HAVE_LANL_INTERPOLATION_LIBRARY
	case NCPA::interpolator1d_t::LANL_1D_LINEAR:
		return static_cast<NCPA::Interpolator1D*>(new NCPA::LANLLinearInterpolator1D());
		break;
	case NCPA::interpolator1d_t::LANL_1D_NATURAL_CUBIC:
		return static_cast<NCPA::Interpolator1D*>(new NCPA::LANLNaturalCubicSplineInterpolator1D());
		break;
#endif

#ifdef HAVE_GSL_INTERPOLATION_LIBRARY
	case NCPA::interpolator1d_t::GSL_1D_LINEAR:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_linear, t));
		break;
	case NCPA::interpolator1d_t::GSL_1D_POLYNOMIAL:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_polynomial, t));
		break;
	case NCPA::interpolator1d_t::GSL_1D_CUBIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_cspline, t));
		break;
	case NCPA::interpolator1d_t::GSL_1D_CUBIC_PERIODIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLPeriodicInterpolator1D(gsl_interp_cspline_periodic, t));
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_akima, t));
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA_PERIODIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLPeriodicInterpolator1D(gsl_interp_akima_periodic, t));
		break;

#if GSL_MAJOR_VERSION > 1
	case NCPA::interpolator1d_t::GSL_1D_STEFFEN:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_steffen, t));
		break;
#endif

#endif

	default:
		throw std::out_of_range("Unrecognized or unimplemented interpolator requested.");
	}
}

NCPA::Interpolator1D::Interpolator1D() : type_{NCPA::interpolator1d_t::INVALID} {}

NCPA::Interpolator1D::Interpolator1D( NCPA::interpolator1d_t t ) : type_{t} {}

NCPA::Interpolator1D::Interpolator1D( const Interpolator1D &other ) {
	type_ = other.type_;
}

NCPA::Interpolator1D::Interpolator1D( Interpolator1D &&other ) {
	::swap(*this,other);
}

void swap( NCPA::Interpolator1D &a, NCPA::Interpolator1D &b ) {
	using std::swap;
	swap( a.type_, b.type_ );
}

/*
 * public interface of evaluation methods:
 * 1. check limits
 * 2. if outside limits, either extrapolate linearly or throw exception
 * 3. otherwise, call proper evaluation method
 */
void NCPA::Interpolator1D::is_extrapolating( bool tf ) {
	extrapolating_ = tf;

}

bool NCPA::Interpolator1D::is_extrapolating() const {
	return extrapolating_;
}

double NCPA::Interpolator1D::f( double x ) {
	return this->df( 0, x );
//
//	if (x < min_x) {
//		if (is_extrapolating()) {
//			double dx_ = x - min_x;
//			double dfdx_ = extrapolate_df_( min_x );
//			double f_ = extrapolate_f_( min_x );
//			return f_ + dfdx_ * dx_;
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else if (x > max_x) {
//		if (is_extrapolating()) {
//			double dx_ = x - max_x;
//			double dfdx_ = extrapolate_df_( max_x );
//			double f_ = extrapolate_f_( max_x );
//			return f_ + dfdx_ * dx_;
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else {
//		return extrapolate_f_( x );
//	}
}

double NCPA::Interpolator1D::df( double x ) {
	return this->df( 1, x );

//	if (x < min_x) {
//		if (is_extrapolating()) {
//			return extrapolate_df_( min_x );
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else if (x > max_x) {
//		if (is_extrapolating()) {
//			return extrapolate_df_( max_x );
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else {
//		return extrapolate_df_( x );
//	}
}

double NCPA::Interpolator1D::df( size_t n, double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x || x > max_x) {
		return extrapolate_( n, x );
	} else {
		return interpolate_( n, x );
	}
}

std::complex<double> NCPA::Interpolator1D::cf( double x ) {
	return this->cdf( 0, x );
}

std::complex<double> NCPA::Interpolator1D::cdf( double x ) {
	return this->cdf( 1, x );
//	double min_x, max_x;
//	get_interp_limits( min_x, max_x );
//
//	if (x < min_x) {
//		if (is_extrapolating()) {
//			return extrapolate_cdf_( min_x );
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else if (x > max_x) {
//		if (is_extrapolating()) {
//			return extrapolate_cdf_( max_x );
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else {
//		return extrapolate_cdf_( x );
//	}
}

std::complex<double> NCPA::Interpolator1D::cdf( size_t n, double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x || x > max_x) {
		return extrapolate_c_( n, x );
	} else {
		return interpolate_c_( n, x );
	}

//	if (n == 0) {
//		return this->cf( x );
//	} else if (n == 1) {
//		return this->cdf( x );
//	}
//	std::complex<double> c0( 0.0, 0.0 );
//
//	// 2nd and higher derivatives are zero for linear extrapolation
//	double min_x, max_x;
//	get_interp_limits( min_x, max_x );
//
//	if (x < min_x) {
//		if (is_extrapolating()) {
//			return c0;
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else if (x > max_x) {
//		if (is_extrapolating()) {
//			return c0;
//		} else {
//			std::ostringstream oss;
//			oss << "Requested value " << x << " outside range of interpolator ["
//					<< min_x << "," << max_x << "]";
//			throw std::out_of_range(oss.str());
//		}
//	} else {
//		return extrapolate_cdf_( n, x );
//	}
}







double NCPA::Interpolator1D::interpolate_f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::interpolate_df_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::interpolate_d2f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::interpolate_d3f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::interpolate_( size_t n, double x ) {
	switch(n) {
	case 0:
		return this->interpolate_f_( x );
		break;
	case 1:
		return this->interpolate_df_( x );
		break;
	case 2:
		return this->interpolate_d2f_( x );
		break;
	case 3:
		return this->interpolate_d3f_( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [0,3].");
	}
}
std::complex<double> NCPA::Interpolator1D::interpolate_cf_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::interpolate_cdf_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::interpolate_cd2f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::interpolate_cd3f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::interpolate_c_( size_t n, double x ) {
	switch(n) {
	case 0:
		return this->interpolate_cf_( x );
		break;
	case 1:
		return this->interpolate_cdf_( x );
		break;
	case 2:
		return this->interpolate_cd2f_( x );
		break;
	case 3:
		return this->interpolate_cd3f_( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [0,3].");
	}
}

const NCPA::interpolator1d_t NCPA::Interpolator1D::type() const {
	return type_;
}

const std::string NCPA::Interpolator1D::identifier() const {
	return NCPA::Interpolator1D::as_string( this->type() );
}


double NCPA::Interpolator1D::extrapolate_f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::extrapolate_df_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::extrapolate_d2f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::extrapolate_d3f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::extrapolate_( size_t n, double x ) {
	switch(n) {
	case 0:
		return this->extrapolate_f_( x );
		break;
	case 1:
		return this->extrapolate_df_( x );
		break;
	case 2:
		return this->extrapolate_d2f_( x );
		break;
	case 3:
		return this->extrapolate_d3f_( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [0,3].");
	}
}
std::complex<double> NCPA::Interpolator1D::extrapolate_cf_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::extrapolate_cdf_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::extrapolate_cd2f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::extrapolate_cd3f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::extrapolate_c_( size_t n, double x ) {
	switch(n) {
	case 0:
		return this->extrapolate_cf_( x );
		break;
	case 1:
		return this->extrapolate_cdf_( x );
		break;
	case 2:
		return this->extrapolate_cd2f_( x );
		break;
	case 3:
		return this->extrapolate_cd3f_( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [0,3].");
	}
}

double NCPA::Interpolator1D::linear_extrapolate_f_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	if (x < min_x) {
		double dx_ = x - min_x;
		double dfdx_ = interpolate_df_( min_x );
		double f_ = interpolate_f_( min_x );
		return f_ + dfdx_ * dx_;
	} else if (x > max_x) {
		double dx_ = x - max_x;
		double dfdx_ = interpolate_df_( max_x );
		double f_ = interpolate_f_( max_x );
		return f_ + dfdx_ * dx_;
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " is within range of interpolator ["
				<< min_x << "," << max_x << "], but extrapolation requested.  "
				<< "You shouldn't be here, something has gone wrong.";
		throw std::logic_error(oss.str());
	}
}

double NCPA::Interpolator1D::linear_extrapolate_df_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	if (x < min_x) {
		return interpolate_df_( min_x );
	} else if (x > max_x) {
		return interpolate_df_( max_x );
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " is within range of interpolator ["
				<< min_x << "," << max_x << "], but extrapolation requested.  "
				<< "You shouldn't be here, something has gone wrong.";
		throw std::logic_error(oss.str());
	}
}

double NCPA::Interpolator1D::linear_extrapolate_d2f_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return 0.0;
}

double NCPA::Interpolator1D::linear_extrapolate_d3f_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return 0.0;
}

std::complex<double> NCPA::Interpolator1D::linear_extrapolate_cf_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	if (x < min_x) {
		double dx_ = x - min_x;
		std::complex<double> dfdx_ = interpolate_cdf_( min_x );
		std::complex<double> f_ = interpolate_cf_( min_x );
		return f_ + dfdx_ * dx_;
	} else if (x > max_x) {
		double dx_ = x - max_x;
		std::complex<double> dfdx_ = interpolate_cdf_( max_x );
		std::complex<double> f_ = interpolate_cf_( max_x );
		return f_ + dfdx_ * dx_;
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " is within range of interpolator ["
				<< min_x << "," << max_x << "], but extrapolation requested.  "
				<< "You shouldn't be here, something has gone wrong.";
		throw std::logic_error(oss.str());
	}
}

std::complex<double> NCPA::Interpolator1D::linear_extrapolate_cdf_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	if (x < min_x) {
		return interpolate_cdf_( min_x );
	} else if (x > max_x) {
		return interpolate_cdf_( max_x );
	} else {
		std::ostringstream oss;
		oss << "Requested value " << x << " is within range of interpolator ["
				<< min_x << "," << max_x << "], but extrapolation requested.  "
				<< "You shouldn't be here, something has gone wrong.";
		throw std::logic_error(oss.str());
	}
}

std::complex<double> NCPA::Interpolator1D::linear_extrapolate_cd2f_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return std::complex<double>(0.0,0.0);
}

std::complex<double> NCPA::Interpolator1D::linear_extrapolate_cd3f_( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return std::complex<double>(0.0,0.0);
}

NCPA::Interpolator1D* NCPA::Interpolator1D::set(
		std::vector<double> x, std::vector<double> y ) {
	if (x.size() != y.size()) {
		throw std::logic_error("Vectors are not of the same size");
	}
	double *xbuffer = NCPA::zeros<double>( x.size() );
	double *ybuffer = NCPA::zeros<double>( y.size() );
	std::copy( x.cbegin(), x.cend(), xbuffer );
	std::copy( y.cbegin(), y.cend(), ybuffer );
	this->set( x.size(), xbuffer, ybuffer );
	delete [] xbuffer;
	delete [] ybuffer;
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::Interpolator1D::set(
		std::vector<double> x, std::vector<std::complex<double>> y ) {
	if (x.size() != y.size()) {
		throw std::logic_error("Vectors are not of the same size");
	}
	double *xbuffer = NCPA::zeros<double>( x.size() );
	std::complex<double> *ybuffer = NCPA::zeros<std::complex<double>>( y.size() );
	std::copy( x.cbegin(), x.cend(), xbuffer );
	std::copy( y.cbegin(), y.cend(), ybuffer );
	this->set( x.size(), xbuffer, ybuffer );
	delete [] xbuffer;
	delete [] ybuffer;
	return static_cast<NCPA::Interpolator1D *>( this );
}
