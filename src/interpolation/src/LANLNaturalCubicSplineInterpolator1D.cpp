#include "Interpolator1D.h"

#ifdef HAVE_LANL_INTERPOLATION_LIBRARY

#include "LANLNaturalCubicSplineInterpolator1D.h"
#include "LANLInterpolation.h"
#include "NCPACommon.h"
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>

NCPA::LANLNaturalCubicSplineInterpolator1D::LANLNaturalCubicSplineInterpolator1D()
		: NCPA::Interpolator1D(), ready_{false}, minx_{0.0}, maxx_{0.0} {
	this->init();
}

NCPA::LANLNaturalCubicSplineInterpolator1D::LANLNaturalCubicSplineInterpolator1D( const LANLNaturalCubicSplineInterpolator1D &other )
		: NCPA::Interpolator1D( other ), ready_{false}, minx_{0.0}, maxx_{0.0} {

	this->init();
	LANL::prep( real_spline_, other.real_spline_.length );
	LANL::prep( imag_spline_, other.imag_spline_.length );

	if (other.real_spline_.x_vals != nullptr) {
		std::copy( other.real_spline_.x_vals,
				other.real_spline_.x_vals +  other.real_spline_.length,
				this->real_spline_.x_vals );
	}
	if (other.real_spline_.f_vals != nullptr) {
		std::copy( other.real_spline_.f_vals,
				other.real_spline_.f_vals +  other.real_spline_.length,
				this->real_spline_.f_vals );
	}
	if (other.imag_spline_.x_vals != nullptr) {
		std::copy( other.imag_spline_.x_vals,
				other.imag_spline_.x_vals +  other.imag_spline_.length,
				this->imag_spline_.x_vals );
	}
	if (other.imag_spline_.f_vals != nullptr) {
		std::copy( other.imag_spline_.f_vals,
				other.imag_spline_.f_vals +  other.imag_spline_.length,
				this->imag_spline_.f_vals );
	}
	minx_ = other.minx_;
	maxx_ = other.maxx_;
	ready();
}

NCPA::LANLNaturalCubicSplineInterpolator1D::LANLNaturalCubicSplineInterpolator1D( LANLNaturalCubicSplineInterpolator1D &&other )
	: NCPA::Interpolator1D() {
	::swap(*this,other);
}

NCPA::LANLNaturalCubicSplineInterpolator1D::~LANLNaturalCubicSplineInterpolator1D() {
	this->free();
}

void swap( NCPA::LANLNaturalCubicSplineInterpolator1D &a, NCPA::LANLNaturalCubicSplineInterpolator1D &b ) {
	using std::swap;
	::swap( static_cast<NCPA::Interpolator1D&>(a),
			static_cast<NCPA::Interpolator1D&>(b) );
	swap( a.real_spline_, b.real_spline_ );
	swap( a.imag_spline_, b.imag_spline_ );
	swap( a.minx_, b.minx_ );
	swap( a.maxx_, b.maxx_ );
	swap( a.ready_, b.ready_ );
}

NCPA::Interpolator1D* NCPA::LANLNaturalCubicSplineInterpolator1D::init() {
	init_spline_( real_spline_ );
	init_spline_( imag_spline_ );
	minx_ = 0;
	maxx_ = 0;
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::init_spline_(LANL::Spline1DNatural &spline) {
	spline.length = 0;
	spline.accel = -1;
	spline.x_vals = nullptr;
	spline.f_vals = nullptr;
	spline.slopes = nullptr;
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::LANLNaturalCubicSplineInterpolator1D::allocate( size_t n ) {
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::allocate_spline_(
		LANL::Spline1DNatural &spline, size_t n ) {
	LANL::prep( spline, n );
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::LANLNaturalCubicSplineInterpolator1D::ready() {
	LANL::set( real_spline_ );
	LANL::set( imag_spline_ );

	if (real_spline_.length > 0) {
		ready_ = true;
	}
	return static_cast<NCPA::Interpolator1D *>( this );
}

bool NCPA::LANLNaturalCubicSplineInterpolator1D::is_ready() {
	return ready_;
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::free() {
	free_spline_( real_spline_ );
	free_spline_( imag_spline_ );
	minx_ = 0.0;
	maxx_ = 0.0;
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::free_spline_(LANL::Spline1DNatural &spline) {
	if (spline.x_vals != nullptr) {
		delete [] spline.x_vals;
		spline.x_vals = nullptr;
	}
	if (spline.f_vals != nullptr) {
		delete [] spline.f_vals;
		spline.f_vals = nullptr;
	}
	if (spline.slopes != nullptr) {
		delete [] spline.slopes;
		spline.slopes = nullptr;
	}
	spline.length = 0;
	spline.accel = 0;
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::LANLNaturalCubicSplineInterpolator1D::set( size_t n, const double *x,
		const double *y ) {
	this->free();
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );

	this->set_x( n, x );
	this->set_y( n, y );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::LANLNaturalCubicSplineInterpolator1D::set( size_t n, const double *x,
		const double *y_real, const double *y_imag ) {
	this->free();
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );

	this->set_x( n, x );
	this->set_y( n, y_real, y_imag );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::LANLNaturalCubicSplineInterpolator1D::set( size_t n, const double *x,
		const std::complex<double> *y ) {
	this->free();
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );

	this->set_x( n, x );
	this->set_y( n, y );
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::set_x( size_t n, const double *x ) {
	if ((size_t)(real_spline_.length) != n) {
		std::ostringstream oss;
		oss << "Spline is set for arrays of length " << real_spline_.length
				<< " but was passed an array of length " << n << ".";
		throw std::out_of_range(oss.str());
	}
	std::copy( x, x+n, real_spline_.x_vals );
	std::copy( x, x+n, imag_spline_.x_vals );
	minx_ = x[0];
	maxx_ = x[n-1];
	ready_ = false;
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::set_y( size_t n, const double *y ) {
	if ((size_t)(real_spline_.length) != n) {
		std::ostringstream oss;
		oss << "Spline is set for arrays of length " << real_spline_.length
				<< " but was passed an array of length " << n << ".";
		throw std::out_of_range(oss.str());
	}
	std::copy( y, y+n, real_spline_.f_vals );
	std::fill( imag_spline_.f_vals, imag_spline_.f_vals+n, 0.0 );
	ready_ = false;
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::set_y( size_t n, const double *y_real,
		const double *y_imag) {
	if ((size_t)(real_spline_.length) != n) {
		std::ostringstream oss;
		oss << "Spline is set for arrays of length " << real_spline_.length
				<< " but was passed an array of length " << n << ".";
		throw std::out_of_range(oss.str());
	}
	std::copy( y_real, y_real+n, real_spline_.f_vals );
	std::copy( y_imag, y_imag+n, imag_spline_.f_vals );
	ready_ = false;
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::set_y( size_t n, const std::complex<double> *y ) {
	if ((size_t)(real_spline_.length) != n) {
		std::ostringstream oss;
		oss << "Spline is set for arrays of length " << real_spline_.length
				<< " but was passed an array of length " << n << ".";
		throw std::out_of_range(oss.str());
	}
	for (size_t i = 0; i < n; i++) {
		real_spline_.f_vals[i] = y[i].real();
		imag_spline_.f_vals[i] = y[i].imag();
	}
	ready_ = false;
}

const std::string NCPA::LANLNaturalCubicSplineInterpolator1D::identifier() const {
	return "LANL 1-D Natural Cubic Spline Interpolator";
}

void NCPA::LANLNaturalCubicSplineInterpolator1D::get_interp_limits( double &xmin, double &xmax ) const {
	xmin = get_low_interp_limit();
	xmax = get_high_interp_limit();
}

double NCPA::LANLNaturalCubicSplineInterpolator1D::get_low_interp_limit() const {
	if (this->real_spline_.x_vals == nullptr || this->imag_spline_.x_vals == nullptr) {
		throw std::out_of_range("Limits for interpolator have not been set");
	}
	return minx_;
}
double NCPA::LANLNaturalCubicSplineInterpolator1D::get_high_interp_limit() const {
	if (this->real_spline_.x_vals == nullptr || this->imag_spline_.x_vals == nullptr) {
		throw std::out_of_range("Limits for interpolator have not been set");
	}
	return maxx_;
}

LANL::Spline1DNatural *NCPA::LANLNaturalCubicSplineInterpolator1D::get_real_spline() {
	return &real_spline_;
}

LANL::Spline1DNatural *NCPA::LANLNaturalCubicSplineInterpolator1D::get_imag_spline() {
	return &imag_spline_;
}

size_t NCPA::LANLNaturalCubicSplineInterpolator1D::max_derivative() const {
	return 3;
}


double NCPA::LANLNaturalCubicSplineInterpolator1D::eval_f_( double x ) {
	if (this->is_ready()) {
		return LANL::eval_f( x, real_spline_ );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
double NCPA::LANLNaturalCubicSplineInterpolator1D::eval_df_( double x ) {
	if (this->is_ready()) {
		return LANL::eval_df( x, real_spline_ );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
double NCPA::LANLNaturalCubicSplineInterpolator1D::eval_d2f_( double x ) {
	if (this->is_ready()) {
		return LANL::eval_ddf( x, real_spline_ );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
double NCPA::LANLNaturalCubicSplineInterpolator1D::eval_d3f_( double x ) {
	if (this->is_ready()) {
		return LANL::eval_dddf( x, real_spline_ );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
std::complex<double> NCPA::LANLNaturalCubicSplineInterpolator1D::eval_cf_( double x ) {
	if (this->is_ready()) {
		double impart = 0.0;
		if (imag_spline_.f_vals != nullptr) {
			impart = LANL::eval_f( x, imag_spline_ );
		}
		return std::complex<double>(
				LANL::eval_f( x, real_spline_ ),
				impart
		);
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
std::complex<double> NCPA::LANLNaturalCubicSplineInterpolator1D::eval_cdf_( double x ) {
	if (this->is_ready()) {
		double impart = 0.0;
		if (imag_spline_.f_vals != nullptr) {
			impart = LANL::eval_df( x, imag_spline_ );
		}
		return std::complex<double>(
				LANL::eval_df( x, real_spline_ ),
				impart
		);
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
std::complex<double> NCPA::LANLNaturalCubicSplineInterpolator1D::eval_cd2f_( double x ) {
	if (this->is_ready()) {
		double impart = 0.0;
		if (imag_spline_.f_vals != nullptr) {
			impart = LANL::eval_ddf( x, imag_spline_ );
		}
		return std::complex<double>(
				LANL::eval_ddf( x, real_spline_ ),
				impart
		);
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
std::complex<double> NCPA::LANLNaturalCubicSplineInterpolator1D::eval_cd3f_( double x ) {
	if (this->is_ready()) {
		double impart = 0.0;
		if (imag_spline_.f_vals != nullptr) {
			impart = LANL::eval_dddf( x, imag_spline_ );
		}
		return std::complex<double>(
				LANL::eval_dddf( x, real_spline_ ),
				impart
		);
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}

#endif