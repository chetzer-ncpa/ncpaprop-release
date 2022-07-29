#include "Turbulence.h"
#include <cassert>
#include <cmath>
#include <complex>
#include "util.h"
#include <random>
#include <vector>
#include <stdexcept>

NCPA::Turbulence::Turbulence( size_t N ) {
	N_ = N;
	G_.reserve( N_ );
	kvec_.reserve( N_ );
	alpha_.reserve( N_ );
	phases_.reserve( N_ );
	set_defaults();
}

NCPA::Turbulence::~Turbulence() { }

void NCPA::Turbulence::set_defaults() {
	// set_reference_temperature( 293.0 );
	set_turbulence_scale( 100.0 );
	set_temperature_factor( 1.0e-10 );
	set_velocity_factor( 1.0e-8 );
}

size_t NCPA::Turbulence::size() const {
	return N_;
}

void NCPA::Turbulence::set_wavenumbers( const std::vector<double> &k_in ) {
	if (k_in.size() != kvec_.size()) {
		throw std::runtime_error( "Size mismatch in Turbulence wavenumbers" );
	}

	kvec_ = k_in;
}

void NCPA::Turbulence::set_wavenumbers_linear( double k1, double k2 ) {
	assert( k1 != k2 );
	if (k2 < k1) {
		double temp = k1;
		k1 = k2;
		k2 = temp;
	}

	double *ls = NCPA::zeros<double>( N_ );
	NCPA::linspace( k1, k2, N_, ls );
	kvec_.assign( ls, ls + N_ );
	delete [] ls;
}

void NCPA::Turbulence::set_wavenumbers_log( double k1, double k2 ) {
	assert( k1 != k2 );
	if (k2 < k1) {
		double temp = k1;
		k1 = k2;
		k2 = temp;
	}

	double *ls = NCPA::logspace<double>( k1, k2, N_ );
	kvec_.assign( ls, ls + N_ );
	delete [] ls;
}

// supplied phases in the range [0,2*PI)
void NCPA::Turbulence::compute_phases( const std::vector<double> &phases ) {
	if (kvec_.size() != N_) {
		throw std::runtime_error( "Wavenumber vector not set for turbulence phase calculation" );
	}

	for (size_t i = 0; i < N_; i++) {
		std::complex<double> ktemp( kvec_[ i ], 0.0 );
		std::complex<double> ptemp( 0.0, phases[ i ] );
		phases_[ i ] = ktemp * std::exp( ptemp );
	}
}

void NCPA::Turbulence::compute_phases() {
	if (kvec_.size() != N_) {
		throw std::runtime_error( "Wavenumber vector not set for turbulence phase calculation" );
	}
	std::vector<double> randn = NCPA::random_numbers( N_, 2.0*PI );
	compute_phases( randn );
}

void NCPA::Turbulence::set_alpha( const std::vector<double> &avec ) {
	if (avec.size() != N_) {
		throw std::runtime_error( "Size mismatch in alpha assignment for turbulence" );
	}
	alpha_ = avec;
}

void NCPA::Turbulence::set_alpha() {
	alpha_ = random_numbers( N_, 2.0*PI );
}

// void NCPA::Turbulence::set_reference_temperature( double T ) {
// 	T0_ = T;
// }

void NCPA::Turbulence::set_turbulence_scale( double Lt ) {
	Lt_ = Lt;
}

void NCPA::Turbulence::set_temperature_factor( double t ) {
	C_T_factor_ = t;
}

void NCPA::Turbulence::set_velocity_factor( double t ) {
	C_v_factor_ = t;
}

void NCPA::Turbulence::compute() {
	// precompute G(k_n)=sqrt(4*pi*deltak*F(k_n)*k_n)
	// from Salomons, p. 226, as used in J.24 and J.27
	std::vector<double> F = von_karman_spectral_density();
	G_[ 0 ] = std::sqrt( 4.0 * PI * F[ 0 ] * kvec_[ 0 ] );
	for (size_t i = 1; i < N_; i++) {
		double delk = kvec_[i] - kvec_[i-1];
		G_[ i ] = std::sqrt( 4.0 * PI * F[ i ] * delk * kvec_[ i ] );
	}
}

double NCPA::Turbulence::get_G( size_t i ) const {
	return G_[ i ];
}

std::complex<double> NCPA::Turbulence::get_k( size_t i ) const {
	return phases_[ i ];
}

double NCPA::Turbulence::get_alpha( size_t i ) const {
	return alpha_[ i ];
}

// Salomons, Eq. I.54
std::vector<double> NCPA::Turbulence::von_karman_spectral_density() {

	double A = 3.300539063636137e-02;
	double B = 1.682618526390545e+00;
	double D = 4.588959617428759e-01;
	double E = 1.223722564647669e+00;
	double G = 1.833333333333333e+00;
	// double c_0 = 331.5*std::sqrt(1 + (T0_ - 273.15)/273.15);
	double K_0 = (2.0*PI)/Lt_;
	// double C_T = std::sqrt( T0_ * T0_ * C_T_factor_ );
	// double C_v = std::sqrt( T0_ * T0_ * C_v_factor_  );
	double CToT02 = C_T_factor_;
	double Cvoc02 = C_v_factor_;
	double ksq, factor1, factor2, factor3;
	std::vector<double> F( N_ );

	for (size_t i = 0; i < N_; i++) {
		ksq = std::pow( phases_[i].real(), 2.0 )
			+ std::pow( phases_[i].imag(), 2.0 );
		factor1 = A / std::pow( (ksq + K_0*K_0), 4.0/3.0 );
		// factor2 = B * std::pow( C_T / T0_, 2.0 ) / 4.0;
		factor2 = B * CToT02 / 4.0;
		factor3 = E * std::pow( phases_[i].imag(), 2.0 ) / (ksq + K_0*K_0);
		// F[ i ] = factor1 * (factor2 + (D + factor3)*G*std::pow(C_v/c_0,2.0));
		F[ i ] = factor1 * (factor2 + (D + factor3)*G*Cvoc02);
	}
	return F;
}