#include "NCPALinearInterpolator1D.h"

#include "Interpolator1D.h"
#include "NCPACommon.h"
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>

NCPA::NCPALinearInterpolator1D::NCPALinearInterpolator1D()
		: NCPA::Interpolator1D(NCPA::interpolator1d_t::NCPA_1D_LINEAR),
		  n_{0} {
	this->init();
}

NCPA::NCPALinearInterpolator1D::NCPALinearInterpolator1D(
		const NCPA::NCPALinearInterpolator1D &other )
		: NCPA::Interpolator1D( other ) {

	this->init();
	if (other.n_ > 0) {
		this->allocate( other.n_ );
	}
	if (other.x_ != nullptr) {
		std::copy(other.x_, other.x_ + other.n_, this->x_ );
	}
	if (other.yr_ != nullptr) {
		std::copy(other.yr_, other.yr_ + other.n_, this->yr_ );
	}
	if (other.yi_ != nullptr) {
		std::copy(other.yi_, other.yi_ + other.n_, this->yi_ );
	}
	if (other.last_interval_ != nullptr) {
		std::copy(other.last_interval_, other.last_interval_ + 2, this->last_interval_ );
	}
}

NCPA::NCPALinearInterpolator1D::NCPALinearInterpolator1D(
		NCPA::NCPALinearInterpolator1D &&other )
	: NCPA::Interpolator1D() {
	::swap(*this,other);
}

NCPA::NCPALinearInterpolator1D::~NCPALinearInterpolator1D() {
	this->free();
}

void swap( NCPA::NCPALinearInterpolator1D &a, NCPA::NCPALinearInterpolator1D &b ) {
	using std::swap;
	::swap( static_cast<NCPA::Interpolator1D&>(a),
			static_cast<NCPA::Interpolator1D&>(b) );
	swap( a.n_, b.n_ );
	swap( a.x_, b.x_ );
	swap( a.yr_, b.yr_ );
	swap( a.yi_, b.yi_ );
	swap( a.last_interval_, b.last_interval_ );
}

NCPA::NCPALinearInterpolator1D& NCPA::NCPALinearInterpolator1D::operator=(
		NCPA::NCPALinearInterpolator1D other ) {
	::swap(*this,other);
	return *this;
}

NCPA::Interpolator1D* NCPA::NCPALinearInterpolator1D::set(
		size_t n, const double *x, const double *y ) {
	if (n != n_ || this->x_ == nullptr || this->yr_ == nullptr || this->yi_ == nullptr) {
		this->init()->allocate(n);
	}
	std::copy( x, x+n, this->x_ );
	std::copy( y, y+n, this->yr_ );

	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPALinearInterpolator1D::set(
		size_t n, const double *x, const double *y_real, const double *y_imag ) {
	if (n != n_ || this->x_ == nullptr || this->yr_ == nullptr || this->yi_ == nullptr) {
		this->init()->allocate(n);
	}
	std::copy( x, x+n, this->x_ );
	std::copy( y_real, y_real+n, this->yr_ );
	std::copy( y_imag, y_imag+n, this->yi_ );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPALinearInterpolator1D::set(
		size_t n, const double *x, const std::complex<double> *y ) {
	if (n != n_ || this->x_ == nullptr || this->yr_ == nullptr || this->yi_ == nullptr) {
		this->init()->allocate(n);
	}
	std::copy( x, x+n, this->x_ );
	complex2double( n, y, this->yr_, this->yi_ );
	return static_cast<NCPA::Interpolator1D *>( this );
}



NCPA::Interpolator1D* NCPA::NCPALinearInterpolator1D::init() {
	free();
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPALinearInterpolator1D::allocate( size_t n ) {
	n_ = n;
	allocate_array_( this->x_, n );
	allocate_array_( this->yr_, n );
	allocate_array_( this->yi_, n );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPALinearInterpolator1D::ready() {
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::NCPALinearInterpolator1D::free() {
	n_ = 0;
	reset_array_( this->x_ );
	reset_array_( this->yr_ );
	reset_array_( this->yi_ );
}

bool NCPA::NCPALinearInterpolator1D::is_ready() const {
	return ((n_ > 0) && (x_ != nullptr) && (yr_ != nullptr) && (yi_ != nullptr)
			&& (x_[n_-1] != x_[0]));
}

size_t NCPA::NCPALinearInterpolator1D::max_derivative() const {
	return 3;
}

void NCPA::NCPALinearInterpolator1D::get_interp_limits( double &min, double &max ) const {
	min = get_low_interp_limit();
	max = get_high_interp_limit();
}

double NCPA::NCPALinearInterpolator1D::get_low_interp_limit() const {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	return x_[0];
}

double NCPA::NCPALinearInterpolator1D::get_high_interp_limit() const {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	return x_[ n_ - 1 ];
}

bool NCPA::NCPALinearInterpolator1D::get_interval_( double x, size_t &i0, size_t &i1 ) {
	bool use_last_interval;
	if (last_interval_ == nullptr) {
		last_interval_ = new double[2];
		use_last_interval = false;
	} else if (x >= last_interval_[0] && x <= last_interval_[1]) {
		i0 = last_interval_[0];
		i1 = last_interval_[1];
		use_last_interval = true;
	} else {
		use_last_interval = false;
	}
	if (use_last_interval) {
		i0 = last_interval_[0];
		i1 = last_interval_[1];
		return true;
	} else {
		return NCPA::find_interval_inclusive<double>( this->x_, this->n_, x, i0, i1);
	}
}

double NCPA::NCPALinearInterpolator1D::interpolate_f_( double x ) {
	size_t i0, i1;
	if (this->get_interval_( x, i0, i1 )) {
		return NCPA::linearInterp( this->x_[i0], this->yr_[i0],
								   this->x_[i1], this->yr_[i1],
								   x );
	} else {
		std::ostringstream oss;
		oss << "Value " << x << " not in interpolatable range [" << x_[0] << ","
				<< x_[n_-1] << "], and extrapolation is turned off." << std::endl;
		throw std::out_of_range( oss.str() );
	}
}

std::complex<double> NCPA::NCPALinearInterpolator1D::interpolate_cf_( double x ) {
	size_t i0, i1;
	if (this->get_interval_( x, i0, i1 )) {
		return std::complex<double>(
			NCPA::linearInterp( this->x_[i0], this->yr_[i0],
								this->x_[i1], this->yr_[i1],
								x ),
			NCPA::linearInterp( this->x_[i0], this->yi_[i0],
								this->x_[i1], this->yi_[i1],
								x )
		);
	} else {
		std::ostringstream oss;
		oss << "Value " << x << " not in interpolatable range [" << x_[0] << ","
				<< x_[n_-1] << "], and extrapolation is turned off." << std::endl;
		throw std::out_of_range( oss.str() );
	}
}

double NCPA::NCPALinearInterpolator1D::interpolate_df_( double x ) {
	size_t i0, i1;
	if (this->get_interval_( x, i0, i1 )) {
		return (this->yr_[i1] - this->yr_[i0]) / (this->x_[i1] - this->x_[i0]);
	} else {
		std::ostringstream oss;
		oss << "Value " << x << " not in interpolatable range [" << x_[0] << ","
				<< x_[n_-1] << "], and extrapolation is turned off." << std::endl;
		throw std::out_of_range( oss.str() );
	}
}

double NCPA::NCPALinearInterpolator1D::interpolate_d2f_( double x ) {
	return 0.0;
}

double NCPA::NCPALinearInterpolator1D::interpolate_d3f_( double x ) {
	return 0.0;
}

std::complex<double> NCPA::NCPALinearInterpolator1D::interpolate_cdf_( double x ) {
	size_t i0, i1;
	if (this->get_interval_( x, i0, i1 )) {
		return std::complex<double>(
			(this->yr_[i1] - this->yr_[i0]) / (this->x_[i1] - this->x_[i0]),
			(this->yi_[i1] - this->yi_[i0]) / (this->x_[i1] - this->x_[i0])
		);
	} else {
		std::ostringstream oss;
		oss << "Value " << x << " not in interpolatable range [" << x_[0] << ","
				<< x_[n_-1] << "], and extrapolation is turned off." << std::endl;
		throw std::out_of_range( oss.str() );
	}
}

std::complex<double> NCPA::NCPALinearInterpolator1D::interpolate_cd2f_( double x ) {
	return std::complex<double>( 0.0, 0.0 );
}

std::complex<double> NCPA::NCPALinearInterpolator1D::interpolate_cd3f_( double x ) {
	return std::complex<double>( 0.0, 0.0 );
}

void NCPA::NCPALinearInterpolator1D::reset_array_( double *&x ) {
	if (x != nullptr) {
		delete [] x;
		x = nullptr;
	}
}

void NCPA::NCPALinearInterpolator1D::allocate_array_( double *&x, size_t n ) {
	x = NCPA::zeros<double>( n );
}

double NCPA::NCPALinearInterpolator1D::extrapolate_f_( double x ) {
	return this->linear_extrapolate_f_( x );
}
double NCPA::NCPALinearInterpolator1D::extrapolate_df_( double x ) {
	return this->linear_extrapolate_df_( x );
}
double NCPA::NCPALinearInterpolator1D::extrapolate_d2f_( double x ) {
	return this->linear_extrapolate_d2f_( x );
}
double NCPA::NCPALinearInterpolator1D::extrapolate_d3f_( double x ) {
	return this->linear_extrapolate_d3f_( x );
}

std::complex<double> NCPA::NCPALinearInterpolator1D::extrapolate_cf_( double x ) {
	return this->linear_extrapolate_cf_( x );
}
std::complex<double> NCPA::NCPALinearInterpolator1D::extrapolate_cdf_( double x ) {
	return this->linear_extrapolate_cdf_( x );
}
std::complex<double> NCPA::NCPALinearInterpolator1D::extrapolate_cd2f_( double x ) {
	return this->linear_extrapolate_cd2f_( x );
}
std::complex<double> NCPA::NCPALinearInterpolator1D::extrapolate_cd3f_( double x ) {
	return this->linear_extrapolate_cd3f_( x );
}

