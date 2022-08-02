#ifndef NCPAPROP_WMODESOLVER_H_INCLUDED
#define NCPAPROP_WMODESOLVER_H_INCLUDED

#ifndef NCPAPROP_WMOD_FILENAME_1D
#define NCPAPROP_WMOD_FILENAME_1D "tloss_1d.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_1D_LOSSLESS
#define NCPAPROP_WMOD_FILENAME_1D_LOSSLESS "tloss_1d.lossless.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_2D
#define NCPAPROP_WMOD_FILENAME_2D "tloss_2d.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_2D_LOSSLESS
#define NCPAPROP_WMOD_FILENAME_2D_LOSSLESS "tloss_2d.lossless.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_1D_MULTIPROP
#define NCPAPROP_WMOD_FILENAME_1D_MULTIPROP "Nby2D_tloss_1d.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_1D_LOSSLESS_MULTIPROP
#define NCPAPROP_WMOD_FILENAME_1D_LOSSLESS_MULTIPROP "Nby2D_tloss_1d.lossless.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_PHASE_SPEEDS
#define NCPAPROP_WMOD_FILENAME_PHASE_SPEEDS "phasespeeds.wnm"
#endif

#ifndef NCPAPROP_WMOD_FILENAME_ATMOSPHERE
#define NCPAPROP_WMOD_FILENAME_ATMOSPHERE "atm_profile.wnm"
#endif


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
