#ifndef NCPAPROP_MODESOLVER_H_INCLUDED
#define NCPAPROP_MODESOLVER_H_INCLUDED

#include "AtmosphericTransferFunctionSolver.h"
#include "parameterset.h"
#include "Atmosphere1D.h"
#include <string>
#include <complex>

// @todo Create destructor and handle memory deallocation there
namespace NCPA {
	class ModeSolver : public AtmosphericTransferFunctionSolver {

	public:           

        // constructor/destructor
		//ModeSolver( ParameterSet *param, Atmosphere1D *atm_profile );
		virtual ~ModeSolver();

		// parameter control
		virtual void setParams( ParameterSet *param, Atmosphere1D *atm_prof ) = 0;                  	
		virtual void printParams() = 0;
		virtual void calculateAttenuation( double freq );

		// general functions
		int getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev);
		int sturmCount(int n, double dz, double *diag, double k, int *cnt);	
		int doPerturb(int nz, double z_min, double dz, int n_modes, double freq,
			double *k, double **v, double *alpha, std::complex<double> *k_pert);
		int getTLoss1D(int select_modes, double dz, int n_r, double dr, double z_src, double z_rcv, 
			double *rho, std::complex<double> *k_pert, double **v_s, std::string filename_lossy, 
			std::string filename_lossless = "");
		int getTLoss1DNx2(double azimuth, int select_modes, double dz, int n_r, double dr, double z_src, 
			double z_rcv,  double *rho, std::complex<double> *k_pert, double **v_s, bool Nx2, int iter,
			std::string filename_lossy, std::string filename_lossless = "" );
		int getTLoss2D(int nz, int select_modes, double dz, int n_r, double dr, double z_src, 
			double *rho, std::complex<double> *k_pert, double **v_s, std::string filename_lossy,
			std::string filename_lossless = "" );
		int writeDispersion(int select_modes, double dz, double z_src, double z_rcv, double freq, 
			std::complex<double> *k_pert, double **v_s, double *rho);
		int writeDispersion(FILE *fp, int select_modes, double dz, double z_src, double z_rcv, 
			double freq, std::complex<double> *k_pert, double **v_s, double *rho);
		int writePhaseSpeeds(int select_modes, double freq,
			std::complex<double> *k_pert, const std::string &filename );
		int writeEigenFunctions(int nz, int select_modes, double dz,
			double **v_s, const std::string &file_extension );

		std::string tag_filename( std::string basename );

		// to be overridden
		virtual int doSelect( int nz, int n_modes, double k_min, double k_max, double *k2, double **v, 
			double *k_s, double **v_s, int *select_modes ) = 0;


	protected:
		bool   write_2D_TLoss;
		bool   write_phase_speeds;
		//bool   write_dispersion;
		bool   write_modes;        
		bool   Nby2Dprop;
		bool   turnoff_WKB;
		bool   wvnum_filter_flg;
		bool   z_min_specified;
		bool   append_dispersion_file;
		bool   write_atmosphere;
		bool   write_lossless;
          
		int    Nz_grid;
		int    Nrng_steps;
		int    Lamb_wave_BC; 
		int    Naz;
		int    Nfreq;
		int    skiplines;
      
		//double freq;
		double azi;
		double azi_min;
		double azi_max;
		double azi_step;           
		double z_min;
		double maxheight;     
		double maxrange;
		double sourceheight;
		double receiverheight;
		double tol;
		double *Hgt, *zw, *mw, *T, *rho, *Pr, *c_eff, *alpha;    // change AA to alpha
		double c_min; // for wavenumber filtering option
		double c_max; // for wavenumber filtering option

		double *f_vec;
//      	double freq;
      	std::string user_tag = "";
      	
		NCPA::Atmosphere1D *atm_profile;
		std::string gnd_imp_model;
		std::string atmosfile;
		//std::string wind_units;
		std::string usrattfile;
		std::string dispersion_file;
      
	}; 
}

#endif
