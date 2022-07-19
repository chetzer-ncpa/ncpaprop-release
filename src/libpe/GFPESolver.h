#ifndef NCPAPROP_GFPESOLVER_H_INCLUDED
#define NCPAPROP_GFPESOLVER_H_INCLUDED

#include <fstream>
#include <complex>
#include <vector>
#include <map>
#include <utility>
#include "util.h"
#include "fftw3.h"
#include "matrix.h"
#include "ncpaprop_common.h"
#include "ncpaprop_atmosphere.h"
#include "Turbulence.h"

// default file names
#define NCPAPROP_GFPE_TLOSS_FILENAME_1D "tloss_1d.pe"
#define NCPAPROP_GFPE_TLOSS_FILENAME_2D "tloss_2d.pe"

namespace NCPA {

	class TransmissionLossField {

	public:
		TransmissionLossField(
			double az, double freq,
			size_t nr,
			double *r_vec,
			size_t nz,
			double *z_vec,
			std::complex<double> **TL
			);

		TransmissionLossField(
			double az, double freq,
			const std::vector<double> &r_vec,
			const std::vector<double> &z_vec,
			const NCPA::DenseMatrix<std::complex<double>> *TL_mat
			);

		~TransmissionLossField();

		// void as_arrays(
		// 	double &az,
		// 	double &freq,
		// 	size_t &nr,
		// 	double *&r_vec,
		// 	size_t &nz,
		// 	double *&z_vec,
		// 	std::complex<double> **&TL
		// 	);

		void as_vectors(
			double &az,
			double &freq,
			std::vector<double> &r_vec,
			std::vector<double> &z_vec,
			NCPA::DenseMatrix<std::complex<double>> *&TL_mat
		) const;

	protected:
		double azimuth_, frequency_;
		std::vector<double> r_vec_, z_vec_;
		NCPA::DenseMatrix<std::complex<double>> *TL_mat_;
	};

	class GFPESolver : public AtmosphericTransferFunctionSolver {

	public:
		GFPESolver();
		GFPESolver( NCPA::ParameterSet *param );
		virtual ~GFPESolver();
		virtual int solve();
		virtual size_t runs() const;
		virtual void truncateFile( const std::string &filename ) const;
		virtual void truncate1DFile() const;
		virtual void truncate2DFile() const;

		virtual void output1DTL(
			const std::string &filename,
			bool append = false
			) const;

		virtual void output1DTL(
			double height,
			const std::string &filename,
			bool append = false
			) const;

		// use default filename
		virtual void output1DTL(
			bool append = false
			) const;

		virtual void output1DTL(
			double height,
			bool append = false
			) const;

		virtual void get1DTL(
			size_t n,
			double height,
			std::vector<double> &r,
			std::vector<std::complex<double>> &tl
			) const;


		virtual void get1DTL(
			size_t n,
			std::vector<double> &r,
			std::vector<std::complex<double>> &tl
			) const;

		virtual void get1DTL(
			double height,
			std::vector<double> &r,
			std::vector<std::complex<double>> &tl
			) const;

		virtual void get1DTL(
			std::vector<double> &r,
			std::vector<std::complex<double>> &tl
			) const;

		virtual void output2DTL(
			const std::string &filename,
			bool append = false
			) const;

		virtual void output2DTL(
			bool append = false
			) const;

		virtual void output2DTL(
			size_t n,
			bool append = false
			) const;

		virtual void output2DTL(
			size_t n,
			const std::string &filename,
			bool append = false
			) const;

		virtual void get2DTL(
			size_t n,
			std::vector<double> &r,
			std::vector<double> &z,
			NCPA::DenseMatrix<std::complex<double>> *&tl_mat ) const;

		virtual void get2DTL(
			std::vector<double> &r,
			std::vector<double> &z,
			NCPA::DenseMatrix<std::complex<double>> *&tl_mat ) const;

	protected:

		void set_default_values();
		void calculate_effective_sound_speed(
			double az, const std::string &new_key );
		double average_value( const std::string &key,
			double range, double low, double high ) const;
		double boundary_layer(
			double z, double A, double ztop, double zatt ) const;
		void wavenumber_filter( size_t nz, double *k, double kzero,
			double f1, double f2, double *&F ) const;
		std::complex<double> ground_impedence( double f );
		void phase_factor( size_t nz, double deltax,
			std::complex<double> *vark, std::complex<double> *&pf );
		void apply_turbulence( size_t nz,
			double k_a, double r, double *z, double *&Gamma ) const;
		void compute_starter( size_t Ntrans, double deltak,
			double k0, double *kprime, std::complex<double> *E3, double *F,
			double *z, double z0, std::complex<double> Zhat, double deltar,
			std::complex<double> *&phi ) const;
		std::string tag_filename( const std::string &base ) const;

		NCPA::Atmosphere2D *atm;           // atmospheric model
		std::vector<double> f_vec, az_vec; // frequencies & azimuths
		bool multiprop;                    // use multiple azimuths?
		bool use_atm_1d, use_atm_2d;       // atmosphere flags
		bool nodelay;					   // include propagation delay
		NCPA::ScalarWithUnits *range;      // propagation distance
		NCPA::ScalarWithUnits *z_source;   // source height
		NCPA::ScalarWithUnits *z_receiver; // receiver height
		NCPA::ScalarWithUnits *z_ground;   // ground height
		NCPA::Turbulence *turbulence;      // turbulence object
		double sigma;						// ground impedence factor
		double deltaz;						// vertical resolution
		int nz_requested;
		double dr_requested;
		double boundary_thickness_requested;
		double surface_thickness_requested;
		double f1, f2;
		std::string filetag;

		// turbulence parameters
		bool use_turbulence, random_turbulence;
		double turbulence_k1, turbulence_k2, Lt,
			temperature_factor, velocity_factor;
		size_t turbulence_size;
		std::string turbulence_file;

		double ablthicknessfactor;			// ABL thickness factor
		double zatmfactor;
		double deltazfactor;
		double deltarfactor;

		bool verbose;

		// store calculation results
		std::vector<NCPA::TransmissionLossField *> TL_;

		// for suite (large-N) runs
		bool use_suite_mode;
		size_t N_suite;
		std::vector<double> r_suite, z_suite;

	};
}






#endif