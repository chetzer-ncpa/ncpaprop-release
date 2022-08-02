#ifndef NCPAPROP_ESSMODESOLVER_H_INCLUDED
#define NCPAPROP_ESSMODESOLVER_H_INCLUDED

#ifndef NCPAPROP_MODESS_FILENAME_1D
#define NCPAPROP_MODESS_FILENAME_1D "tloss_1d.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_1D_LOSSLESS
#define NCPAPROP_MODESS_FILENAME_1D_LOSSLESS "tloss_1d.lossless.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_2D
#define NCPAPROP_MODESS_FILENAME_2D "tloss_2d.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_2D_LOSSLESS
#define NCPAPROP_MODESS_FILENAME_2D_LOSSLESS "tloss_2d.lossless.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_1D_MULTIPROP
#define NCPAPROP_MODESS_FILENAME_1D_MULTIPROP "Nby2D_tloss_1d.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_1D_LOSSLESS_MULTIPROP
#define NCPAPROP_MODESS_FILENAME_1D_LOSSLESS_MULTIPROP "Nby2D_tloss_1d.lossless.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_PHASE_SPEEDS
#define NCPAPROP_MODESS_FILENAME_PHASE_SPEEDS "phasespeeds.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_PHASE_AND_GROUP_SPEEDS
#define NCPAPROP_MODESS_FILENAME_PHASE_AND_GROUP_SPEEDS "speeds.nm"
#endif

#ifndef NCPAPROP_MODESS_FILENAME_ATMOSPHERE
#define NCPAPROP_MODESS_FILENAME_ATMOSPHERE "atm_profile.nm"
#endif




#include "slepceps.h"
#include "slepcst.h"
#include "ModeSolver.h"
#include "parameterset.h"
#include "Atmosphere1D.h"
#include <string>
#include <complex>

namespace NCPA {

	class ESSModeSolver : public ModeSolver {

	public:
		ESSModeSolver( ParameterSet *param, Atmosphere1D *atm_profile );
		~ESSModeSolver();

		void setParams( ParameterSet *param, Atmosphere1D *atm_prof );
		void printParams();

		int solve();
		//int doESSSLEPcCalculation( double *diag, double dz, double *k_min, double *k_max, 
		//	PetscInt *nconv, double *k2, double **v );
		int doSelect( int nz, int n_modes, double k_min, double k_max, double *k2, double **v, 
			double *k_s, double **v_s, int *select_modes );
		int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, 
			double dz, NCPA::Atmosphere1D *p, double admittance, double freq, double azi, 
			double *diag, double *k_min, double *k_max, bool turnoff_WKB, double *c_eff);
		void getModalStarter(int nz, int select_modes, double dz, double freq,  double z_src, 
			double z_rcv, double *rho, std::complex<double> *k_pert, double **v_s, 
			const std::string &modstartfile);
		int writePhaseAndGroupSpeeds(int nz, double dz, int select_modes, double freq, 
			std::complex<double> *k_pert, double **v_s, double *c_eff,
			const std::string &filename );

	protected:
		std::string modstartfile; // store the modal starter in this file
		bool   write_speeds;
	};
}




#endif
