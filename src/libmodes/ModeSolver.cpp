#include <complex>
#include <stdexcept>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "AtmosphericTransferFunctionSolver.h"
#include "modes.h"
#include "Atmosphere1D.h"
#include "util.h"

#define MAX_MODES 4000 

using namespace NCPA;
using namespace std;

NCPA::ModeSolver::~ModeSolver() {
	delete[] Hgt;
	delete[] zw;
	delete[] mw;
	delete[] T;
	delete[] rho;
	delete[] Pr;
	delete[] c_eff;
	delete[] alpha;
	delete[] f_vec;

}

void NCPA::ModeSolver::calculateAttenuation( double freq ) {
	atm_profile->remove_property( "_ALPHA_" );
	if (usrattfile.empty()) {
		atm_profile->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", freq );
	} else {
		atm_profile->read_attenuation_from_file( "_ALPHA_", usrattfile );
	}
	
	for (int i=0; i<Nz_grid; i++) {
		alpha[i]   = atm_profile->get( "_ALPHA_", Hgt[i] );
	}
}



int NCPA::ModeSolver::getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev)
{
	int nev_max, nev_min;
	sturmCount(n,dz,diag,k_max,&nev_max);
	sturmCount(n,dz,diag,k_min,&nev_min);
	*nev = nev_max - nev_min;
	return 0;
}


int NCPA::ModeSolver::sturmCount(int n, double dz, double *diag, double k, int *cnt)
{
	double kk,pot,cup0,cup1,cup2;
	double fd_d_val, fd_o_val;
	int i,pm;

	fd_d_val = -2./pow(dz,2);  // Finite Difference Coefficient on  diagonal
	fd_o_val =  1./pow(dz,2);  // Finite Difference Coefficient off diagonal

	pm   = 0;
	kk   = k*k;
	cup0 = fd_d_val + diag[n-1] - kk;
	pot  = fd_d_val + diag[n-2] - kk;
	cup1 = cup0*pot;
	if (cup0*cup1 < 0.0) { pm++; } 
	cup0=cup0/fabs(cup1);
	cup1=cup1/fabs(cup1);

	for (i=n-3; i>=0; i--) {
		pot  = fd_d_val + diag[i] - kk;
		cup2 = pot*cup1 - (pow(fd_o_val,2))*cup0;
		if (cup1*cup2 < 0.0) { pm++; }
		cup0=cup1/fabs(cup2);
		cup1=cup2/fabs(cup2);
	}
	*cnt = pm;
	return 0;
}


int NCPA::ModeSolver::doPerturb(int nz, double z_min, double dz, int n_modes, double freq, /*SampledProfile *p,*/ 
		double *k, double **v, double *alpha, complex<double> *k_pert) {
	int i, j;
	double absorption, gamma, c_T;
	double z_km, dz_km;
	double omega = 2*PI*freq;
	complex<double> I (0.0, 1.0);    // @todo replace with I macro from <complex>
	gamma = 1.4;
	
	dz_km = dz/1000.0;
	for (j=0; j<n_modes; j++) {
		absorption = 0.0;
		z_km=z_min/1000.0;
		for (i=0; i<nz; i++) {
			c_T = sqrt(gamma*Pr[i]/rho[i]); // in m/s  @todo create and get from c0 vector?
			absorption = absorption + dz*v[i][j]*v[i][j]*(omega/c_T)*alpha[i]*2;
			z_km += dz_km;
		}			
		k_pert[j] = sqrt(k[j]*k[j] + I*absorption);
	}
	return 0;
}






// 20150603 DV added code to account for sqrt(rho_rcv/rho_src) if so chosen
int NCPA::ModeSolver::getTLoss1D(int select_modes, double dz, int n_r, double dr, double z_src, 
		double z_rcv, double *rho, complex<double> *k_pert, double **v_s, 
		std::string filename_lossy, std::string filename_lossless ) {

	// computes the transmissison loss from source to receiver.
	// The formula for pressure is:
	// p(r,z) = sqrt(rho(z)/rho(zs))*I*exp(-I*pi/4)/sqrt(8*pi*r)*sum(Vm(z)*Vm(zs)*exp(I*k_m*r)/sqrt(k_m))
	// where Vm are the modes computed in this code. Note that Vm(z) = Psi_m(z)/sqrt(rho(z)), 
	// where Psi_m(z) are the modes defined in Ocean Acoustics, 1994 ed. pg. 274.
	// Note that we usually save the TL corresponding to
	// the reduced pressure formula: p_red(r,z) = p(r,z)/sqrt(rho(z))
	// check below to see which formula is actually implemented.
			
	int i, m;
	int n_zsrc = (int) ceil(z_src/dz);
	int n_zrcv = (int) ceil(z_rcv/dz);
	double r, sqrtrho_ratio;
	double modal_sum_i, modal_sum_i_ll; // for incoherent sum if desired
	complex<double> modal_sum_c, modal_sum_c_ll;
	complex<double> I (0.0, 1.0);
	
	// the 4*PI factor ensures that the modal sum below ends up being the actual TL
	complex<double> expov8pi =  4*PI*I*exp(-I*PI*0.25)/sqrt(8.0*PI); 
  
	sqrtrho_ratio  = sqrt(rho[n_zrcv]/rho[n_zsrc]);
  
    FILE *tloss_1d, *tloss_ll_1d;
	tloss_1d    = fopen( filename_lossy.c_str(), "w");
	if (write_lossless) {
		tloss_ll_1d = fopen( filename_lossless.c_str(), "w");
	}
	// @todo check that files were opened properly

	for (i=0; i<n_r; i++) {
		r = (i+1)*dr;
		modal_sum_c    = 0.;
		modal_sum_c_ll = 0;
		modal_sum_i    = 0;      
		modal_sum_i_ll = 0;
      
		// here we save the modal sum; the reduced pressure is:
		//     p_red(r,z) =modal_sum/(4*pi*sqrt(rho(zs)))= p(r,z)/sqrt(rho(z))
		for (m=0; m<select_modes; m++) {    
			modal_sum_c    = modal_sum_c + v_s[n_zsrc][m] * v_s[n_zrcv][m]
					 * exp(I*k_pert[m]*r) / sqrt(k_pert[m]);
			modal_sum_i    = modal_sum_i + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)
					 * exp(-2*imag(k_pert[m])*r) / abs(k_pert[m]);

			if (write_lossless) {
				modal_sum_c_ll = modal_sum_c_ll + v_s[n_zsrc][m] * v_s[n_zrcv][m]
					 * exp(I*real(k_pert[m])*r) / sqrt(real(k_pert[m]));
				modal_sum_i_ll = modal_sum_i_ll + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)
					 / real(k_pert[m]);
			}

		}
      
		// no sqrt(rho[n_zrcv]/rho[n_zsrc]) factor
		// if (1) {
		modal_sum_c    = expov8pi*modal_sum_c/sqrt(r);
		modal_sum_i    = 4*PI*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
		if (write_lossless) {
			modal_sum_c_ll = expov8pi*modal_sum_c_ll/sqrt(r);
			modal_sum_i_ll = 4*PI*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		}
		// }
      
		// sqrt(rho[n_zrcv]/rho[n_zsrc]) factor added
		// if (0) {
		// 	modal_sum_c    = sqrtrho_ratio*expov8pi*modal_sum_c/sqrt(r);
		// 	modal_sum_c_ll = sqrtrho_ratio*expov8pi*modal_sum_c_ll/sqrt(r);
		// 	modal_sum_i    = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
		// 	modal_sum_i_ll = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		// }

		fprintf(tloss_1d,"%f %20.12e %20.12e %20.12e\n", r/1000.0, real(modal_sum_c), 
			imag(modal_sum_c), modal_sum_i);
		if (write_lossless) {
			fprintf(tloss_ll_1d,"%f %20.12e %20.12e %20.12e\n", r/1000.0, real(modal_sum_c_ll),
				imag(modal_sum_c_ll), modal_sum_i_ll);
		}
	}
	fclose(tloss_1d);
	printf("           file %s created\n", filename_lossy.c_str() );
	if (write_lossless) {
		fclose(tloss_ll_1d);
		printf("           file %s created\n", filename_lossless.c_str() );
	}
	return 0;
}




// 20150603 DV added code to account for sqrt(rho_rcv/rho_src) if so chosen
int NCPA::ModeSolver::getTLoss1DNx2(double azimuth, int select_modes, double dz, int n_r, double dr, 
	double z_src, double z_rcv, double *rho, complex<double> *k_pert, double **v_s, 
	bool Nx2, int iter, std::string filename_lossy, std::string filename_lossless ) {

	// !!! note the modes assumed here are the modes in Oc. Acoust. divided by sqrt(rho(z))
	// Thus the formula for (physical) pressure is:
	// p(r,z) = i*exp(-i*pi/4)/sqrt(8*pi*r)*1/rho(zs)* sum[ v_s(zs)*sqrt(rho(zs))*v_s(z)*sqrt(rho(z))
	// *exp(ikr)/sqrt(k) ]
	// We usually save the TL corresponding to the reduced pressure: p_red(r,z) = p(r,z)/sqrt(rho(z))

	int i, m;
	int n_zsrc = (int) ceil(z_src/dz);
	int n_zrcv = (int) ceil(z_rcv/dz);
	double r, sqrtrho_ratio;
	double modal_sum_i, modal_sum_i_ll;
	complex<double> modal_sum_c, modal_sum_c_ll;
	complex<double> I (0.0, 1.0);
	
	// the 4*PI factor ensures that the modal sum below ends up being the actual TL
	complex<double> expov8pi =  4*PI*I*exp(-I*PI*0.25)/sqrt(8.0*PI); 
	FILE *tloss_1d, *tloss_ll_1d;
  
	sqrtrho_ratio = sqrt(rho[n_zrcv]/rho[n_zsrc]);
	
	if (iter==0) {
		tloss_1d    = fopen(filename_lossy.c_str(),"w");
		if (write_lossless) {
			tloss_ll_1d = fopen(filename_lossless.c_str(),"w");
		}
	}
	else {  // append
		tloss_1d    = fopen(filename_lossy.c_str(),"a");
		if (write_lossless) {
			tloss_ll_1d = fopen(filename_lossless.c_str(),"a");
		}
	}
	// @todo make sure files were opened properly

	for (i=0; i<n_r; i++) {
		r = (i+1)*dr;
		modal_sum_c = 0.;
		modal_sum_i = 0.;
		modal_sum_c_ll = 0;
		modal_sum_i_ll = 0;
      
		// Use the commented lines if the modes must be scaled by sqrt(rho)
		// @todo under what circumstances should we activate these?
		/*
		for (m=0; m<select_modes; m++) {
		modal_sum_c    = modal_sum_c + v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio*exp(I*k_pert[m]*r)/sqrt(k_pert[m]);
		modal_sum_i    = modal_sum_i + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio,2)*exp(-2*imag(k_pert[m])*r)/abs(k_pert[m]);
		modal_sum_c_ll = modal_sum_c_ll + v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio*exp(I*real(k_pert[m])*r)/sqrt(real(k_pert[m]));
		modal_sum_i_ll = modal_sum_i_ll + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio,2)/real(k_pert[m]);
		}
		*/
      
		// here we save the modal sum; the reduced pressure is: p_red(r,z) =modal_sum/(4*pi*sqrt(rho(zs)))= p(r,z)/sqrt(rho(z))
		for (m=0; m<select_modes; m++) {
			modal_sum_c    = modal_sum_c 
				+ (v_s[n_zsrc][m] * v_s[n_zrcv][m] * exp(I*k_pert[m]*r) / sqrt(k_pert[m]));
			modal_sum_i    = modal_sum_i
				+ pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)*exp(-2*imag(k_pert[m])*r)/abs(k_pert[m]);
			if (write_lossless) {
				modal_sum_c_ll = modal_sum_c_ll
					+ v_s[n_zsrc][m]*v_s[n_zrcv][m]*exp(I*real(k_pert[m])*r)/sqrt(real(k_pert[m]));
				modal_sum_i_ll = modal_sum_i_ll + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)/real(k_pert[m]);
			}
		}
      
      
		// no sqrt(rho[n_zrcv]/rho[n_zsrc]) factor
		// @todo Need programmatic logic to determine which branch to take, or remove one
		modal_sum_c    = expov8pi*modal_sum_c/sqrt(r);
		modal_sum_i    = 4*PI*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
		if (write_lossless) {
			modal_sum_c_ll = expov8pi*modal_sum_c_ll/sqrt(r);
			modal_sum_i_ll = 4*PI*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		}

      
		// sqrt(rho[n_zrcv]/rho[n_zsrc]) factor added
		// if (0) {
		// 	modal_sum_c    = sqrtrho_ratio*expov8pi*modal_sum_c/sqrt(r);
		// 	modal_sum_c_ll = sqrtrho_ratio*expov8pi*modal_sum_c_ll/sqrt(r);
		// 	modal_sum_i    = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
		// 	modal_sum_i_ll = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		// }
      
		fprintf(tloss_1d,"%10.3f %8.3f %20.12e %20.12e %20.12e\n", r/1000.0, 
			azimuth, real(modal_sum_c), imag(modal_sum_c), modal_sum_i);
		if (write_lossless) {
			fprintf(tloss_ll_1d,"%10.3f %8.3f %20.12e %20.12e %20.12e\n", r/1000.0,
				azimuth, real(modal_sum_c_ll), imag(modal_sum_c_ll), modal_sum_i_ll);
		}
	}
	fprintf(tloss_1d, "\n");
	fclose(tloss_1d);
	if (write_lossless) {
		fprintf(tloss_ll_1d, "\n");
		fclose(tloss_ll_1d);
	}
	return 0;
}

// 20150603 DV added code to account for sqrt(rho_rcv/rho_src) if so chosen
int NCPA::ModeSolver::getTLoss2D(int nz, int select_modes, double dz, int n_r, double dr, 
	double z_src, double *rho, complex<double> *k_pert, double **v_s, 
	std::string filename_lossy, std::string filename_lossless ) {

	// !!! note the modes assumed here are the modes in Oc. Acoust. divided by sqrt(rho(z))
	// Thus the formula for (physical) pressure is:
	// p(r,z) = i*exp(-i*pi/4)/sqrt(8*pi*r)*1/rho(zs)* sum[ v_s(zs)*sqrt(rho(zs))*v_s(z)*sqrt(rho(z))
	// *exp(ikr)/sqrt(k) ]
	// We usually save the TL corresponding to the reduced pressure: p_red(r,z) = p(r,z)/sqrt(rho(z))
  
	int i, j, m, stepj;
	int n_zsrc = (int) ceil(z_src/dz);
	double r, z, sqrtrhoj, rho_atzsrc;
	complex<double> modal_sum, modal_sum_ll;
	complex<double> I (0.0, 1.0);
	
	// the 4*PI factor ensures that the modal sum below ends up being the actual TL
	complex<double> expov8pi = 4*PI*I*exp(-I*PI*0.25)/sqrt(8.0*PI); 
  
	rho_atzsrc = rho[n_zsrc];

	stepj = nz/500; // controls the vertical sampling of 2D data saved 
	if (stepj==0) {
		stepj = 10;	// ensure it's never 0; necessary for the loop below
	}

	FILE *tloss_2d, *tloss_2d_ll;
	tloss_2d = fopen(filename_lossy.c_str(),"w");
	if (write_lossless) {
		tloss_2d_ll = fopen(filename_lossless.c_str(), "w" );
	}

	for (i=0; i<n_r; i++) {
		r = (i+1)*dr;
		for (j=0; j<nz; j=j+stepj) {
			z = (j)*dz;
			sqrtrhoj = sqrt(rho[j]);
			modal_sum = 0.0;
			modal_sum_ll = 0.0;
          
			// Use the commented lines if the modes must be scaled by sqrt(rho)
			for (m=0; m<select_modes; m++) {
				modal_sum = modal_sum + v_s[n_zsrc][m]*v_s[j][m]
					*exp(I*k_pert[m]*r)/sqrt(k_pert[m]);
				if (write_lossless) {
					modal_sum_ll += v_s[n_zsrc][m]*v_s[j][m]
						* std::exp(I*k_pert[m].real()*r) / std::sqrt(k_pert[m].real());
				}
			}
			modal_sum = expov8pi*modal_sum/sqrt(r); // no sqrt(rho[n_zrcv]/rho[n_zsrc]) factor
			if (write_lossless) {
				modal_sum_ll *= expov8pi / std::sqrt(r);
			}

			fprintf(tloss_2d,"%f %f %15.8e %15.8e\n", r/1000.0, z/1000.0,
				real(modal_sum), imag(modal_sum));
			if (write_lossless) {
				fprintf(tloss_2d_ll,"%f %f %15.8e %15.8e\n", r/1000.0, z/1000.0,
					real(modal_sum_ll), imag(modal_sum_ll));
			}
		}
		fprintf(tloss_2d,"\n");
		if (write_lossless) {
			fprintf(tloss_2d_ll,"\n");
		}
	}
	fclose(tloss_2d);
	printf("           file %s created\n", filename_lossy.c_str() );
	if (write_lossless) {
		fclose(tloss_2d_ll);
		printf("           file %s created\n", filename_lossless.c_str() );
	}
	return 0;
}


// DV 20150409 - added rho as argument
int NCPA::ModeSolver::writeDispersion(FILE *dispersion, int select_modes, double dz, double z_src, 
	double z_rcv, double freq, complex<double> *k_pert, double **v_s, double *rho) {
	int i;
	int n_zsrc = (int) ceil(z_src/dz);
	int n_zrcv = (int) ceil(z_rcv/dz);
	//char dispersion_file[128];
  
	//printf("In writeDispersion(): n_zsrc=%d ; n_zrcv=%d; z_src=%g; z_rcv=%g\n", n_zsrc, n_zrcv, z_src, z_rcv);

	//sprintf(dispersion_file,"dispersion_%e.nm", freq);
	//FILE *dispersion = fopen(dispersion_file,"w");
	fprintf(dispersion,"%.7e   %d    %.7e   %.7e", freq, select_modes, rho[n_zsrc], rho[n_zrcv]);
	for( i=0; i<select_modes; i++ ) {
		fprintf(dispersion,"   %.12e   %.12e",real(k_pert[i]),imag(k_pert[i]));
		fprintf(dispersion,"   %.12e   %.12e",(v_s[n_zsrc][i]),(v_s[n_zrcv][i]));
	}
	fprintf(dispersion,"\n");
	//fclose(dispersion);
	//printf("           file %s created\n", dispersion_file.c_str());
	return 0;
}

int NCPA::ModeSolver::writeDispersion(int select_modes, double dz, double z_src, 
	double z_rcv, double freq, complex<double> *k_pert, double **v_s, double *rho) {

	char dispersion_file[128];
	sprintf(dispersion_file,"dispersion_%e.nm", freq);
	FILE *dispersion = fopen(dispersion_file,"w");
	int rval = writeDispersion( dispersion, select_modes, dz, z_src, z_rcv, freq, k_pert, v_s, rho );
	fclose( dispersion );
	return rval;
}
 

int NCPA::ModeSolver::writePhaseSpeeds(int select_modes, double freq,
	complex<double> *k_pert, const std::string &filename ) {
	int j;
	FILE *phasespeeds= fopen(filename.c_str(), "w");
	// @todo check to see if file is opened
	for (j=0; j< select_modes; j++) {
		fprintf(phasespeeds, "%d %f %15.8e\n", j, (2*PI*freq)/real(k_pert[j]), imag(k_pert[j]));
	}
	fclose(phasespeeds);
	printf("           file %s created\n",filename.c_str());
	return 0;
}



int NCPA::ModeSolver::writeEigenFunctions(int nz, int select_modes,
	double dz, double **v_s, const std::string &file_extension)
{
	char mode_output[40];
	int j, n;
	double chk;
	double dz_km = dz/1000.0;

	for (j=0; j<select_modes; j++) {
		sprintf(mode_output,"mode_%d.%s", j, file_extension.c_str());
		FILE *eigenfunction= fopen(tag_filename(mode_output).c_str(), "w");
		chk = 0.0;
		for (n=0; n<nz; n++) {
			fprintf(eigenfunction,"%f %15.8e\n", n*dz_km, v_s[n][j]);
			chk = chk + v_s[n][j]*v_s[n][j]*dz;
		}
		if (fabs(1.-chk) > 0.1) { 
			printf("Check if eigenfunction %d is normalized!\n", j); 
		}
		fclose(eigenfunction);
	}
	printf("           %d mode files created\n", select_modes);
	return 0;
}

std::string NCPA::ModeSolver::tag_filename( std::string basename ) {
	return user_tag + basename;
}
