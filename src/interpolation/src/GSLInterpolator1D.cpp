#include "Interpolator1D.h"


// this will only compile if gsl/gsl_spline.h is available
#ifdef HAVE_GSL_INTERPOLATION_LIBRARY

#include "gsl/gsl_spline.h"
#include "GSLInterpolator1D.h"
#include "NCPACommon.h"
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>

NCPA::GSLInterpolator1D::GSLInterpolator1D() : NCPA::Interpolator1D(),
	ready_{false}, minx_{0.0}, maxx_{0.0}, interptype_{nullptr} {}

NCPA::GSLInterpolator1D::GSLInterpolator1D(
		const gsl_interp_type *interp_type,
		NCPA::interpolator1d_t ncpa_type )
		: NCPA::Interpolator1D(ncpa_type), ready_{false}, minx_{0.0}, maxx_{0.0} {
	interptype_ = new gsl_interp_type( *interp_type );
	this->init();
}

NCPA::GSLInterpolator1D::GSLInterpolator1D( const GSLInterpolator1D &other )
		: NCPA::Interpolator1D( other ), ready_{false}, minx_{0.0}, maxx_{0.0} {
	interptype_ = new gsl_interp_type( *(other.get_gsl_interp_type()) );
	this->init();
	if (other.spline_r_ != nullptr && other.spline_i_ != nullptr) {
		allocate_spline_( spline_r_, accel_r_, other.spline_r_->size );
		allocate_spline_( spline_i_, accel_i_, other.spline_i_->size );
		this->set( other.spline_r_->size, other.spline_r_->x, other.spline_r_->y,
				other.spline_i_->y );
		this->accel_r_->cache = other.accel_r_->cache;
		this->accel_i_->cache = other.accel_i_->cache;
		this->accel_r_->hit_count = other.accel_r_->hit_count;
		this->accel_i_->hit_count = other.accel_i_->hit_count;
		this->accel_r_->miss_count = other.accel_r_->miss_count;
		this->accel_i_->miss_count = other.accel_i_->miss_count;
	}
	this->minx_ = other.minx_;
	this->maxx_ = other.maxx_;
	ready();
}

NCPA::GSLInterpolator1D::GSLInterpolator1D( GSLInterpolator1D &&other )
	: NCPA::Interpolator1D() {
	::swap(*this,other);
}

NCPA::GSLInterpolator1D::~GSLInterpolator1D() {
	this->free();
	delete interptype_;
}

void swap( NCPA::GSLInterpolator1D &a, NCPA::GSLInterpolator1D &b ) {
	using std::swap;
	::swap( static_cast<NCPA::Interpolator1D&>(a),
			static_cast<NCPA::Interpolator1D&>(b) );
	swap( a.spline_r_, b.spline_r_ );
	swap( a.spline_i_, b.spline_i_ );
	swap( a.accel_r_, b.accel_r_ );
	swap( a.accel_i_, b.accel_i_ );
	swap( a.interptype_, b.interptype_ );
	swap( a.ready_, b.ready_ );
	swap( a.minx_, b.minx_ );
	swap( a.maxx_, b.maxx_ );
}

NCPA::Interpolator1D* NCPA::GSLInterpolator1D::init() {
	spline_r_ = nullptr;
	spline_i_ = nullptr;
	accel_r_ = nullptr;
	accel_i_ = nullptr;
	minx_ = 0.0;
	maxx_ = 0.0;
	return static_cast<NCPA::Interpolator1D *>( this );
}

//void NCPA::GSLInterpolator1D::init_spline_(gsl_spline * &spline) {
//
//}

NCPA::Interpolator1D* NCPA::GSLInterpolator1D::allocate( size_t n ) {
	allocate_spline_( spline_r_, accel_r_, n );
	allocate_spline_( spline_i_, accel_i_, n );
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::GSLInterpolator1D::allocate_spline_(
		gsl_spline *&spline, gsl_interp_accel *&accel, size_t n ) {
	spline = gsl_spline_alloc( interptype_, n );
	accel  = gsl_interp_accel_alloc();
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::GSLInterpolator1D::ready() {
	if (spline_r_ != nullptr && spline_i_ != nullptr) {
		ready_ = true;
	} else {
		throw std::domain_error("Interpolator has not been set up");
	}
	return static_cast<NCPA::Interpolator1D *>( this );
}

bool NCPA::GSLInterpolator1D::is_ready() const {
	return ready_;
}

void NCPA::GSLInterpolator1D::free() {
	free_spline_( spline_r_, accel_r_ );
	free_spline_( spline_i_, accel_i_ );
	minx_ = 0.0;
	maxx_ = 0.0;
}

void NCPA::GSLInterpolator1D::free_spline_(gsl_spline *&spline,
		gsl_interp_accel *&accel) {
	if (spline != nullptr) {
		gsl_spline_free( spline );
		spline = nullptr;
	}
	if (accel != nullptr) {
		gsl_interp_accel_free( accel );
		accel = nullptr;
	}
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::GSLInterpolator1D::set( size_t n, const double *x, const double *y ) {
	double *z = NCPA::zeros<double>( n );
	this->set_splines_( n, x, y, z );
	delete [] z;
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::GSLInterpolator1D::set( size_t n, const double *x,
		const double *y_r, const double *y_i ) {
//	for (size_t i = 0; i < n; i++) {
//		std::cout << "y_r[" << i << "] = " << y_r[i] << std::endl;
//		std::cout << "y_i[" << i << "] = " << y_i[i] << std::endl;
//	}
	this->set_splines_( n, x, y_r, y_i );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::GSLInterpolator1D::set( size_t n, const double *x,
		const std::complex<double> *y ) {
	double *y_r = new double[n];
	double *y_i = new double[n];
	complex2double( n, y, y_r, y_i );
	this->set_splines_( n, x, y_r, y_i );
	delete [] y_r;
	delete [] y_i;
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::GSLInterpolator1D::set_splines_( size_t n, const double *x, const double *y_r,
		const double *y_i ) {
	if (n < interptype_->min_size) {
		std::ostringstream oss;
		oss << "GSL interpolation type " << interptype_->name
				<< " requires at least " << interptype_->min_size
				<< " points.";
		throw std::domain_error(oss.str());
	}
	this->free();
	allocate_spline_( spline_r_, accel_r_, n );
	allocate_spline_( spline_i_, accel_i_, n );
	gsl_spline_init( spline_r_, x, y_r, n );
	gsl_spline_init( spline_i_, x, y_i, n );
	minx_ = x[0];
	maxx_ = x[n-1];
}

size_t NCPA::GSLInterpolator1D::max_derivative() const {
	return 2;
}


double NCPA::GSLInterpolator1D::interpolate_f_( double x ) {
	if (this->is_ready()) {
		return gsl_spline_eval( spline_r_, x, accel_r_ );
	} else {
		throw std::domain_error( "Interpolator not ready!" );
	}
}
double NCPA::GSLInterpolator1D::interpolate_df_( double x ) {
	if (this->is_ready()) {
		return gsl_spline_eval_deriv( spline_r_, x, accel_r_ );
	} else {
		throw std::domain_error( "Interpolator not ready!" );
	}
}
double NCPA::GSLInterpolator1D::interpolate_d2f_( double x ) {
	if (this->is_ready()) {
		return gsl_spline_eval_deriv2( spline_r_, x, accel_r_ );
	} else {
		throw std::domain_error( "Interpolator not ready!" );
	}
}

std::complex<double> NCPA::GSLInterpolator1D::interpolate_cf_( double x ) {
	if (this->is_ready()) {
		return std::complex<double>(
				gsl_spline_eval( spline_r_, x, accel_r_ ),
				gsl_spline_eval( spline_i_, x, accel_i_ )
		);
	} else {
		throw std::domain_error( "Interpolator not ready!" );
	}
}
std::complex<double> NCPA::GSLInterpolator1D::interpolate_cdf_( double x ) {
	if (this->is_ready()) {
		return std::complex<double>(
				gsl_spline_eval_deriv( spline_r_, x, accel_r_ ),
				gsl_spline_eval_deriv( spline_i_, x, accel_i_ )
		);
	} else {
		throw std::domain_error( "Interpolator not ready!" );
	}
}
std::complex<double> NCPA::GSLInterpolator1D::interpolate_cd2f_( double x ) {
	if (this->is_ready()) {
		return std::complex<double>(
				gsl_spline_eval_deriv2( spline_r_, x, accel_r_ ),
				gsl_spline_eval_deriv2( spline_i_, x, accel_i_ )
		);
	} else {
		throw std::domain_error( "Interpolator not ready!" );
	}
}

gsl_interp_type *NCPA::GSLInterpolator1D::get_gsl_interp_type() const {
	return interptype_;
}

void NCPA::GSLInterpolator1D::set_gsl_interp_type( const gsl_interp_type *gsltype ) {
	interptype_ = new gsl_interp_type(*gsltype);
}

void NCPA::GSLInterpolator1D::get_interp_limits( double &xmin, double &xmax ) const {
	xmin = get_low_interp_limit();
	xmax = get_high_interp_limit();
}

double NCPA::GSLInterpolator1D::get_low_interp_limit() const {
	if (spline_r_ == nullptr) {
		throw std::domain_error( "Interpolator not ready!" );
	}
	return minx_;
}
double NCPA::GSLInterpolator1D::get_high_interp_limit() const {
	if (spline_r_ == nullptr) {
		throw std::domain_error( "Interpolator not ready!" );
	}
	return maxx_;
}

gsl_spline *NCPA::GSLInterpolator1D::get_real_spline() {
	return spline_r_;
}

gsl_spline *NCPA::GSLInterpolator1D::get_imag_spline() {
	return spline_i_;
}

gsl_interp_accel *NCPA::GSLInterpolator1D::get_real_accel() {
	return accel_r_;
}

gsl_interp_accel *NCPA::GSLInterpolator1D::get_imag_accel() {
	return accel_i_;
}

double NCPA::GSLInterpolator1D::extrapolate_f_( double x ) {
	return this->linear_extrapolate_f_( x );
}
double NCPA::GSLInterpolator1D::extrapolate_df_( double x ) {
	return this->linear_extrapolate_df_( x );
}
double NCPA::GSLInterpolator1D::extrapolate_d2f_( double x ) {
	return this->linear_extrapolate_d2f_( x );
}

std::complex<double> NCPA::GSLInterpolator1D::extrapolate_cf_( double x ) {
	return this->linear_extrapolate_cf_( x );
}
std::complex<double> NCPA::GSLInterpolator1D::extrapolate_cdf_( double x ) {
	return this->linear_extrapolate_cdf_( x );
}
std::complex<double> NCPA::GSLInterpolator1D::extrapolate_cd2f_( double x ) {
	return this->linear_extrapolate_cd2f_( x );
}


#endif
