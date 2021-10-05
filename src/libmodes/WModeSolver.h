#ifndef NCPAPROP_WMODESOLVER_H_INCLUDED
#define NCPAPROP_WMODESOLVER_H_INCLUDED

#include "ModeSolver.h"
#include "parameterset.h"
#include "Atmosphere1D.h"
#include <string>
#include <complex>

namespace NCPA {

	class WModeSolver : public ModeSolver {

	public:
		WModeSolver( ParameterSet *param, Atmosphere1D *atm_profile );
		~WModeSolver();

		void setParams( ParameterSet *param, Atmosphere1D *atm_prof );
		void printParams();

		int solve();
		int doSelect( int nz, int n_modes, double k_min, double k_max, double *k2, double **v, 
			double *k_s, double **v_s, int *select_modes );
		int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, double dz, Atmosphere1D *p, 
            double admittance, double freq, double *diag, double *kd, double *md, double *cd, double *k_min, 
            double *k_max, bool turnoff_WKB);

	};
}




#endif
