#ifndef NCPAPROP_EIGENENGINE_H_INCLUDED
#define NCPAPROP_EIGENENGINE_H_INCLUDED

#include "slepceps.h"
#include "slepcst.h"
#include <string>

namespace NCPA {

	class EigenEngine {

	public:
		static int doESSCalculation( double *diag, int Nz_grid, double dz, double tol, 
			double *k_min, double *k_max, PetscInt *nconv, double *k2, double **v );
		static int doWideAngleCalculation( int Nz_grid, double dz, double k_min, double k_max,
			double tol, int nev, double *kd, double *md, double *cd, 
			PetscInt *nconv, double *kH, double **v, std::string disp_msg );


	};

}








#endif