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
#include "GSLPeriodicInterpolator1D.h"
#endif

#include <complex>
#include <stdexcept>
#include <sstream>


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
				new NCPA::GSLPeriodicInterpolator1D(gsl_interp_cspline_periodic));
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLInterpolator1D(gsl_interp_akima));
		break;
	case NCPA::interpolator1d_t::GSL_1D_AKIMA_PERIODIC:
		return static_cast<NCPA::Interpolator1D*>(
				new NCPA::GSLPeriodicInterpolator1D(gsl_interp_akima_periodic));
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

/*
 * public interface of evaluation methods:
 * 1. check limits
 * 2. if outside limits, either extrapolate linearly or throw exception
 * 3. otherwise, call proper evaluation method
 */
void NCPA::Interpolator1D::is_extrapolating( bool tf ) {
	extrapolate_ = tf;

}

bool NCPA::Interpolator1D::is_extrapolating() const {
	return extrapolate_;
}

double NCPA::Interpolator1D::f( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x) {
		if (is_extrapolating()) {
			double dx_ = x - min_x;
			double dfdx_ = eval_df_( min_x );
			double f_ = eval_f_( min_x );
			return f_ + dfdx_ * dx_;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else if (x > max_x) {
		if (is_extrapolating()) {
			double dx_ = x - max_x;
			double dfdx_ = eval_df_( max_x );
			double f_ = eval_f_( max_x );
			return f_ + dfdx_ * dx_;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else {
		return eval_f_( x );
	}
}

double NCPA::Interpolator1D::df( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x) {
		if (is_extrapolating()) {
			return eval_df_( min_x );
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else if (x > max_x) {
		if (is_extrapolating()) {
			return eval_df_( max_x );
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else {
		return eval_df_( x );
	}
}

double NCPA::Interpolator1D::df( size_t n, double x ) {
	if (n == 0) {
		return this->f( x );
	} else if (n == 1) {
		return this->df( x );
	}

	// 2nd and higher derivatives are zero for linear extrapolation
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x) {
		if (is_extrapolating()) {
			return 0.0;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else if (x > max_x) {
		if (is_extrapolating()) {
			return 0.0;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else {
		return eval_df_( n, x );
	}
}

std::complex<double> NCPA::Interpolator1D::cf( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x) {
		if (is_extrapolating()) {
			double dx_ = x - min_x;
			std::complex<double> dfdx_ = eval_cdf_( min_x );
			std::complex<double> f_ = eval_cf_( min_x );
			return f_ + dfdx_ * dx_;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else if (x > max_x) {
		if (is_extrapolating()) {
			double dx_ = x - max_x;
			std::complex<double> dfdx_ = eval_cdf_( max_x );
			std::complex<double> f_ = eval_cf_( max_x );
			return f_ + dfdx_ * dx_;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else {
		return eval_cf_( x );
	}
}

std::complex<double> NCPA::Interpolator1D::cdf( double x ) {
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x) {
		if (is_extrapolating()) {
			return eval_cdf_( min_x );
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else if (x > max_x) {
		if (is_extrapolating()) {
			return eval_cdf_( max_x );
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else {
		return eval_cdf_( x );
	}
}

std::complex<double> NCPA::Interpolator1D::cdf( size_t n, double x ) {
	if (n == 0) {
		return this->cf( x );
	} else if (n == 1) {
		return this->cdf( x );
	}
	std::complex<double> c0( 0.0, 0.0 );

	// 2nd and higher derivatives are zero for linear extrapolation
	double min_x, max_x;
	get_interp_limits( min_x, max_x );

	if (x < min_x) {
		if (is_extrapolating()) {
			return c0;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else if (x > max_x) {
		if (is_extrapolating()) {
			return c0;
		} else {
			std::ostringstream oss;
			oss << "Requested value " << x << " outside range of interpolator ["
					<< min_x << "," << max_x << "]";
			throw std::out_of_range(oss.str());
		}
	} else {
		return eval_cdf_( n, x );
	}
}







double NCPA::Interpolator1D::eval_f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_df_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_d2f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_d3f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
double NCPA::Interpolator1D::eval_df_( size_t n, double x ) {
	switch(n) {
	case 0:
		return this->eval_f_( x );
		break;
	case 1:
		return this->eval_df_( x );
		break;
	case 2:
		return this->eval_d2f_( x );
		break;
	case 3:
		return this->eval_d3f_( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [0,3].");
	}
}
std::complex<double> NCPA::Interpolator1D::eval_cf_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cdf_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cd2f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cd3f_( double x ) {
	throw std::domain_error( "Not implemented for this type of interpolator" );
}
std::complex<double> NCPA::Interpolator1D::eval_cdf_( size_t n, double x ) {
	switch(n) {
	case 0:
		return this->eval_cf_( x );
		break;
	case 1:
		return this->eval_cdf_( x );
		break;
	case 2:
		return this->eval_cd2f_( x );
		break;
	case 3:
		return this->eval_cd3f_( x );
		break;
	default:
		throw std::out_of_range( "Derivative order must be in [0,3].");
	}
}
