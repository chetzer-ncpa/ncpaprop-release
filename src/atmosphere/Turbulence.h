#ifndef NCPA_TURBULENCE_H_INCLUDED
#define NCPA_TURBULENCE_H_INCLUDED

#ifndef PI
#define PI 3.14159265358979323846
#endif

#include <cstdlib>
#include <complex>
#include <vector>

namespace NCPA {

	class Turbulence {

	public:
		Turbulence( size_t N );
		Turbulence( const Turbulence &orig );
		~Turbulence();

		size_t size() const;

		void set_defaults();
		// void set_reference_temperature( double T );
		void set_turbulence_scale( double Lt );
		void set_temperature_factor( double d );
		void set_velocity_factor( double d );

		void compute();

		// set wavenumber spectrum
		void set_wavenumbers(
			const std::vector<double> &kvec );
		void set_wavenumbers_log( double k1, double k2 );
		void set_wavenumbers_linear( double k1, double k2 );

		// compute random phases
		void compute_phases( const std::vector<double> &randnum );
		void compute_phases();

		// do we need these?
		void set_alpha( const std::vector<double> &avec );
		void set_alpha();

		double get_G( size_t i ) const;
		std::complex<double> get_k( size_t ) const;
		double get_alpha( size_t i ) const;

		std::vector<double> von_karman_spectral_density();

	protected:
		size_t N_;
		std::vector<double> G_, kvec_;
		std::vector<std::complex<double>> phases_;
		std::vector<double> alpha_;
		double Lt_, C_T_factor_, C_v_factor_;
	};




	// typedef struct turbulence_struct {
	// 	size_t N;
	// 	double *G;
	// 	std::complex<double> *kvec;
	// 	double *alpha;
	// } turbulence_struct;

	// void init_turbulence_struct( turbulence_struct &ts, size_t N );
	// void free_turbulence_struct( turbulence_struct &ts );

	// void vonKarmanG( size_t n_random_angles,
	// 	double kmin, double kmax, double C_T, double C_v, double T_0,
	// 	double L, turbulence_struct &atmosphere );

	// void logspace( double a, double b, size_t k, double *&ls );

	// void vonKarmanSpectralDensity( size_t N, std::complex<double> *kvec,
	// 	double C_T, double C_v, double T_0, double L, double *&F );

}


#endif
