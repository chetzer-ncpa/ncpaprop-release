#include "NCPANearestNeighborInterpolator1D.h"

#include "Interpolator1D.h"
#include <complex>
#include <algorithm>
#include "NCPACommon.h"



NCPA::NCPANearestNeighborInterpolator1D::NCPANearestNeighborInterpolator1D()
		: NCPA::Interpolator1D(NCPA::interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR) {}

NCPA::NCPANearestNeighborInterpolator1D::NCPANearestNeighborInterpolator1D(
		const NCPA::NCPANearestNeighborInterpolator1D &other ) : NCPA::Interpolator1D(other) {
	this->n_ = other.n_;
	if (other.x_ != nullptr) {
		this->x_ = NCPA::zeros<double>( this->n_ );
		std::copy( other.x_, other.x_ + other.n_, this->x_ );
	}
	if (other.yr_ != nullptr) {
		this->yr_ = NCPA::zeros<double>( this->n_ );
		std::copy( other.yr_, other.yr_ + other.n_, this->yr_ );
	}
	if (other.yi_ != nullptr) {
		this->yi_ = NCPA::zeros<double>( this->n_ );
		std::copy( other.yi_, other.yi_ + other.n_, this->yi_ );
	}
}

NCPA::NCPANearestNeighborInterpolator1D::NCPANearestNeighborInterpolator1D(
		NCPA::NCPANearestNeighborInterpolator1D &&other ) : NCPA::Interpolator1D() {
	::swap( *this, other );
}

void swap( NCPA::NCPANearestNeighborInterpolator1D &a, NCPA::NCPANearestNeighborInterpolator1D &b ) {
	using std::swap;
	::swap(static_cast<NCPA::Interpolator1D &>(a), static_cast<NCPA::Interpolator1D &>(b));
	swap(a.n_, b.n_);
	swap(a.x_, b.x_);
	swap(a.yr_, b.yr_);
	swap(a.yi_, b.yi_);
}

NCPA::NCPANearestNeighborInterpolator1D::~NCPANearestNeighborInterpolator1D() {
	this->init();
}

NCPA::NCPANearestNeighborInterpolator1D& NCPA::NCPANearestNeighborInterpolator1D::operator=(
		NCPA::NCPANearestNeighborInterpolator1D other ) {
	::swap(*this,other);
	return *this;
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::clone() const {
	return static_cast<NCPA::Interpolator1D *>( new NCPA::NCPANearestNeighborInterpolator1D( *this ) );
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::set(
		size_t n, const double *x, const double *y ) {
	if (n != this->n_ || this->x_ == nullptr || this->yr_ == nullptr) {
		this->init()->allocate(n);
	}
	std::copy( x, x+n, this->x_ );
	std::copy( y, y+n, this->yr_ );

	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::set(
		size_t n, const double *x, const double *y_real, const double *y_imag ) {
	if (n != this->n_ || this->x_ == nullptr || this->yr_ == nullptr
			|| this->yi_ == nullptr) {
		this->init()->allocate(n);
	}
	std::copy( x, x+n, this->x_ );
	std::copy( y_real, y_real+n, this->yr_ );
	std::copy( y_imag, y_imag+n, this->yi_ );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::set(
		size_t n, const double *x, const std::complex<double> *y ) {
	if (n != this->n_ || this->x_ == nullptr || this->yr_ == nullptr
			|| this->yi_ == nullptr) {
		this->init()->allocate(n);
	}
	std::copy( x, x+n, this->x_ );
	complex2double( n, y, this->yr_, this->yi_ );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::init() {
	free();
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::allocate( size_t n ) {
	n_ = n;
	allocate_array_( this->x_, n );
	allocate_array_( this->yr_, n );
	allocate_array_( this->yi_, n );
	return static_cast<NCPA::Interpolator1D *>( this );
}

NCPA::Interpolator1D* NCPA::NCPANearestNeighborInterpolator1D::ready() {
	if ((n_ == 0) || (x_ == nullptr) || (yr_ == nullptr) || (x_[n_-1] == x_[0])) {
		throw std::domain_error("Interpolator has not been set up");
	}
	return static_cast<NCPA::Interpolator1D *>( this );
}

void NCPA::NCPANearestNeighborInterpolator1D::free() {
	n_ = 0;
	reset_array_( this->x_ );
	reset_array_( this->yr_ );
	reset_array_( this->yi_ );
}

bool NCPA::NCPANearestNeighborInterpolator1D::is_ready() const {
	return ((n_ > 0) && (x_ != nullptr) && (yr_ != nullptr) && (yi_ != nullptr)
			&& (x_[n_-1] != x_[0]));
}

size_t NCPA::NCPANearestNeighborInterpolator1D::max_derivative() const {
	return 3;
}

void NCPA::NCPANearestNeighborInterpolator1D::get_interp_limits( double &min, double &max ) const {
	min = get_low_interp_limit();
	max = get_high_interp_limit();
}

double NCPA::NCPANearestNeighborInterpolator1D::get_low_interp_limit() const {
	if (!this->is_ready()) {
		throw std::domain_error( "Interpolator not ready!" );
	}
	return x_[0];
}

double NCPA::NCPANearestNeighborInterpolator1D::get_high_interp_limit() const {
	if (!this->is_ready()) {
		throw std::domain_error( "Interpolator not ready!" );
	}
	return x_[ n_ - 1 ];
}

double NCPA::NCPANearestNeighborInterpolator1D::interpolate_f_( double x ) {
	if (this->is_ready()) {
		return yr_[ NCPA::find_closest_index<double>( this->x_, n_, x ) ];
	} else {
		throw std::domain_error("Interpolator not ready");
	}
}

std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::interpolate_cf_( double x ) {
	size_t ind = NCPA::find_closest_index<double>( this->x_, n_, x );
	return std::complex<double>( yr_[ind], yi_[ind] );
}

double NCPA::NCPANearestNeighborInterpolator1D::interpolate_df_( double x ) {
	return 0.0;
}

double NCPA::NCPANearestNeighborInterpolator1D::interpolate_d2f_( double x ) {
	return 0.0;
}

double NCPA::NCPANearestNeighborInterpolator1D::interpolate_d3f_( double x ) {
	return 0.0;
}


std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::interpolate_cdf_( double x ) {
	return std::complex<double>( 0.0, 0.0 );
}

std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::interpolate_cd2f_( double x ) {
	return std::complex<double>( 0.0, 0.0 );
}

std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::interpolate_cd3f_( double x ) {
	return std::complex<double>( 0.0, 0.0 );
}

void NCPA::NCPANearestNeighborInterpolator1D::reset_array_( double *&x ) {
	if (x != nullptr) {
		delete [] x;
		x = nullptr;
	}
}

void NCPA::NCPANearestNeighborInterpolator1D::allocate_array_( double *&x, size_t n ) {
	x = NCPA::zeros<double>( n );
}

//bool NCPA::NCPANearestNeighborInterpolator1D::is_extrapolating() const { return true; }

double NCPA::NCPANearestNeighborInterpolator1D::extrapolate_f_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return yr_[ NCPA::find_closest_index<double>( this->x_, n_, x ) ];
}

std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::extrapolate_cf_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	size_t ind = NCPA::find_closest_index<double>( this->x_, n_, x );
	return std::complex<double>( yr_[ind], yi_[ind] );
}

double NCPA::NCPANearestNeighborInterpolator1D::extrapolate_df_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
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

double NCPA::NCPANearestNeighborInterpolator1D::extrapolate_d2f_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
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

double NCPA::NCPANearestNeighborInterpolator1D::extrapolate_d3f_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
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


std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::extrapolate_cdf_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return std::complex<double>( 0.0, 0.0 );
}

std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::extrapolate_cd2f_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return std::complex<double>( 0.0, 0.0 );
}

std::complex<double> NCPA::NCPANearestNeighborInterpolator1D::extrapolate_cd3f_( double x ) {
	if (!this->is_ready()) {
		throw std::domain_error("Interpolator not ready");
	}
	double min_x, max_x;
	get_interp_limits( min_x, max_x );
	if (!is_extrapolating()) {
		std::ostringstream oss;
		oss << "Requested value " << x << " outside range of interpolator ["
				<< min_x << "," << max_x << "], and extrapolation is turned off.";
		throw std::out_of_range(oss.str());
	}
	return std::complex<double>( 0.0, 0.0 );
}




