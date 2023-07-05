#include "Interpolator1D.h"

#ifdef HAVE_LANL_INTERPOLATION_LIBRARY

#include "LANLLinearInterpolator1D.h"
#include "LANLInterpolation.h"
#include "NCPACommon.h"
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>

NCPA::LANLLinearInterpolator1D::LANLLinearInterpolator1D()
		: NCPA::Interpolator1D(), ready_{false}, minx_{0.0}, maxx_{0.0} {
	this->init();
}

NCPA::LANLLinearInterpolator1D::LANLLinearInterpolator1D( const LANLLinearInterpolator1D &other )
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
	this->minx_ = other.minx_;
	this->maxx_ = other.maxx_;
	ready();
}

NCPA::LANLLinearInterpolator1D::LANLLinearInterpolator1D( LANLLinearInterpolator1D &&other )
	: NCPA::Interpolator1D() {
	::swap(*this,other);
}

NCPA::LANLLinearInterpolator1D::~LANLLinearInterpolator1D() {
	this->free();
}

void swap( NCPA::LANLLinearInterpolator1D &a, NCPA::LANLLinearInterpolator1D &b ) {
	using std::swap;
	::swap( static_cast<NCPA::Interpolator1D&>(a),
			static_cast<NCPA::Interpolator1D&>(b) );
	swap( a.real_spline_, b.real_spline_ );
	swap( a.imag_spline_, b.imag_spline_ );
	swap( a.minx_, b.minx_ );
	swap( a.maxx_, b.maxx_ );
	swap( a.ready_, b.ready_ );
}

NCPA::Interpolator1D* NCPA::LANLLinearInterpolator1D::init() {
	init_spline_( real_spline_ );
	init_spline_( imag_spline_ );
	minx_ = 0.0;
	maxx_ = 0.0;
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::LANLLinearInterpolator1D::init_spline_(LANL::Spline1DLinear &spline) {
	spline.length = 0;
	spline.accel = -1;
	spline.x_vals = nullptr;
	spline.f_vals = nullptr;
	spline.slopes = nullptr;
	minx_ = 0.0;
	maxx_ = 0.0;
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::LANLLinearInterpolator1D::allocate( size_t n ) {
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::LANLLinearInterpolator1D::allocate_spline_(
		LANL::Spline1DLinear &spline, size_t n ) {
	LANL::prep( spline, n );
	ready_ = false;
}

NCPA::Interpolator1D* NCPA::LANLLinearInterpolator1D::ready() {
	LANL::set( real_spline_ );
	LANL::set( imag_spline_ );

	if (real_spline_.length > 0) {
		ready_ = true;
	}
	return static_cast<NCPA::Interpolator1D *>( this );
}

bool NCPA::LANLLinearInterpolator1D::is_ready() {
	return ready_;
}

void NCPA::LANLLinearInterpolator1D::free() {
	free_spline_( real_spline_ );
	free_spline_( imag_spline_ );
	minx_ = 0.0;
	maxx_ = 0.0;
	init();
}

void NCPA::LANLLinearInterpolator1D::free_spline_(LANL::Spline1DLinear &spline) {
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

NCPA::Interpolator1D* NCPA::LANLLinearInterpolator1D::set( size_t n, const double *x, const double *y ) {
	this->free();
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );

	this->set_x( n, x );
	this->set_y( n, y );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::LANLLinearInterpolator1D::set( size_t n, const double *x,
		const double *y_real, const double *y_imag ) {
	this->free();
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );

	this->set_x( n, x );
	this->set_y( n, y_real, y_imag );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::LANLLinearInterpolator1D::set( size_t n, const double *x,
		const std::complex<double> *y ) {
	this->free();
	allocate_spline_( real_spline_, n );
	allocate_spline_( imag_spline_, n );

	this->set_x( n, x );
	this->set_y( n, y );
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::LANLLinearInterpolator1D::set_x( size_t n, const double *x ) {
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

void NCPA::LANLLinearInterpolator1D::set_y( size_t n, const double *y ) {
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

void NCPA::LANLLinearInterpolator1D::set_y( size_t n, const double *y_real,
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

void NCPA::LANLLinearInterpolator1D::set_y( size_t n, const std::complex<double> *y ) {
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

const std::string NCPA::LANLLinearInterpolator1D::identifier() const {
	return "LANL 1-D Linear Interpolator";
}

size_t NCPA::LANLLinearInterpolator1D::max_derivative() const {
	return 3;
}

void NCPA::LANLLinearInterpolator1D::get_interp_limits( double &xmin, double &xmax ) const {
	xmin = get_low_interp_limit();
	xmax = get_high_interp_limit();
}

double NCPA::LANLLinearInterpolator1D::get_low_interp_limit() const {
	if (this->real_spline_.x_vals == nullptr || this->imag_spline_.x_vals == nullptr) {
		throw std::out_of_range("Limits for interpolator have not been set");
	}
	return minx_;
}
double NCPA::LANLLinearInterpolator1D::get_high_interp_limit() const {
	if (this->real_spline_.x_vals == nullptr || this->imag_spline_.x_vals == nullptr) {
		throw std::out_of_range("Limits for interpolator have not been set");
	}
	return maxx_;
}


double NCPA::LANLLinearInterpolator1D::eval_f_( double x ) {
	if (this->is_ready()) {
		return LANL::eval_f( x, real_spline_ );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
double NCPA::LANLLinearInterpolator1D::eval_df_( double x ) {
	if (this->is_ready()) {
		return LANL::eval_df( x, real_spline_ );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
double NCPA::LANLLinearInterpolator1D::eval_d2f_( double x ) {
	if (this->is_ready()) {
		return 0.0;
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
double NCPA::LANLLinearInterpolator1D::eval_d3f_( double x ) {
	if (this->is_ready()) {
		return 0.0;
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}

std::complex<double> NCPA::LANLLinearInterpolator1D::eval_cf_( double x ) {
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
std::complex<double> NCPA::LANLLinearInterpolator1D::eval_cdf_( double x ) {
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
std::complex<double> NCPA::LANLLinearInterpolator1D::eval_cd2f_( double x ) {
	if (this->is_ready()) {
		return std::complex<double>( 0.0, 0.0 );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}
std::complex<double> NCPA::LANLLinearInterpolator1D::eval_cd3f_( double x ) {
	if (this->is_ready()) {
		return std::complex<double>( 0.0, 0.0 );
	} else {
		throw std::domain_error( "Spline not ready!" );
	}
}

LANL::Spline1DLinear * NCPA::LANLLinearInterpolator1D::get_real_spline() {
	return &real_spline_;
}
LANL::Spline1DLinear * NCPA::LANLLinearInterpolator1D::get_imag_spline() {
	return &imag_spline_;
}

#endif
