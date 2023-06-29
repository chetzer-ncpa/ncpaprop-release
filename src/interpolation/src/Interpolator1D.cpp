#include "Interpolator1D.h"

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
#endif

#include <complex>
#include <stdexcept>



NCPA::Interpolator1D *NCPA::Interpolator1D::build( NCPA::interpolator1d_t t ) {
	switch (t) {
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
				new NCPA::GSLInterpolator1D(gsl_interp_linear));
		break;
	case NCPA::interpolator1d_t::GSL_1D_POLYNOMIAL:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_polynomial));
		break;
	case NCPA::interpolator1d_t::GSL_1D_CUBIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_cspline));
		break;
	case NCPA::interpolator1d_t::GSL_1D_CUBIC_PERIODIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_cspline_periodic));
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_akima));
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA_PERIODIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_akima_periodic));
		break;

#if GSL_MAJOR_VERSION > 1
	case NCPA::interpolator1d_t::GSL_1D_STEFFEN:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_steffen));
		break;
#endif

#endif

	default:
		throw std::out_of_range("Unrecognized or unimplemented interpolator requested.");
	}
}


NCPA::Interpolator1D::Interpolator1D() {}

NCPA::Interpolator1D::Interpolator1D( const Interpolator1D &other ) {}

NCPA::Interpolator1D::Interpolator1D( Interpolator1D &&other ) {
	::swap(*this,other);
}

void swap( NCPA::Interpolator1D &a, NCPA::Interpolator1D &b ) {}

double NCPA::Interpolator1D::eval_f( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_df( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_ddf( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_dddf( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_df( size_t n, double x ) {
	switch(n) {
	case 1:
		return this->eval_df( x );
		break;
	case 2:
		return this->eval_ddf( x );
		break;
	case 3:
		return this->eval_dddf( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [1,3].");
	}
}
std::complex<double> NCPA::Interpolator1D::eval_cf( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cdf( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cddf( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cdddf( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cdf( size_t n, double x ) {
	switch(n) {
	case 1:
		return this->eval_cdf( x );
		break;
	case 2:
		return this->eval_cddf( x );
		break;
	case 3:
		return this->eval_cdddf( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [1,3].");
	}
}
