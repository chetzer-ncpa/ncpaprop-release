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

//
// constructor
//
//NCPA::ModeSolver::ModeSolver( NCPA::ParameterSet *param, NCPA::Atmosphere1D *atm_profile )
//{
//	setParams( param, atm_profile );           
//}

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
	// atm_profile->remove_property( "_WS_" );
	// atm_profile->remove_property( "_WD_" );
	// atm_profile->remove_property( "_C0_" );
	// atm_profile->remove_property( "_ALPHA_" );
}

/*
// setParams() prototype
void NCPA::ModeSolver::setParams( NCPA::ParameterSet *param, NCPA::Atmosphere1D *atm_prof )
{		

	// obtain the parameter values from the user's options
	atmosfile 			= param->getString( "atmosfile" );
	gnd_imp_model 		= param->getString( "ground_impedence_model" );
	usrattfile 			= param->getString( "use_attn_file" );
	modstartfile 		= param->getString( "modal_starter_file" );
  	z_min 				= param->getFloat( "zground_km" ) * 1000.0;    // meters
  	freq 				= param->getFloat( "freq" );
  	maxrange 			= param->getFloat( "maxrange_km" ) * 1000.0;
  	maxheight 			= param->getFloat( "maxheight_km" ) * 1000.0;      // @todo fix elsewhere that m is required
  	sourceheight 		= param->getFloat( "sourceheight_km" );
  	receiverheight 		= param->getFloat( "receiverheight_km" );
  	tol 				= 1.0e-8;
  	Nz_grid 			= param->getInteger( "Nz_grid" );
  	Nrng_steps 			= param->getInteger( "Nrng_steps" );
  	Lamb_wave_BC 		= param->getBool( "Lamb_wave_BC" );
  	write_2D_TLoss  	= param->getBool( "write_2D_TLoss" );
  	write_phase_speeds 	= param->getBool( "write_phase_speeds" );
  	write_speeds 		= param->getBool( "write_speeds" );
  	write_modes 		= param->getBool( "write_modes" );
  	write_dispersion 	= param->getBool( "write_dispersion" );
  	Nby2Dprop 			= param->getBool( "Nby2Dprop" );
  	turnoff_WKB 		= param->getBool( "turnoff_WKB" );
  	z_min_specified     = param->wasFound( "zground_km" );

	// default values for c_min, c_max and wvnum_filter_flg
	c_min = 0.0;
	c_max = 0.0;

	// set c_min, c_max if wavenumber filtering is on
	wvnum_filter_flg  	= param->getBool( "wvnum_filter" );  
	if (wvnum_filter_flg==1) {
		c_min = param->getFloat( "c_min" );
		c_max = param->getFloat( "c_max" );
	};


	if (write_phase_speeds || write_speeds || write_2D_TLoss 
			       || write_modes || write_dispersion) {
		turnoff_WKB = 1; // don't use WKB least phase speed estim. when saving any of the above values
	}
  
	Naz = 1; // Number of azimuths: default is a propagation along a single azimuth
	if (Nby2Dprop) {
		azi  		= param->getFloat( "azimuth_start" );
		azi_max 	= param->getFloat( "azimuth_end" );
		azi_step 	= param->getFloat( "azimuth_step" );
		Naz      = (int) ((azi_max - azi)/azi_step) + 1;
	} else {
		azi 				= param->getFloat( "azimuth" );
	}

	azi_min     = azi;
	atm_profile = atm_prof;

	// get Hgt, zw, mw, T, rho, Pr in SI units; they are freed in computeModes()
	// @todo write functions to allocate and free these, it shouldn't be hidden
	// Note units: height in m, wind speeds in m/s, pressure in Pa(?), density in
	// kg/m3
	Hgt   = new double [Nz_grid];
	zw    = new double [Nz_grid];
	mw    = new double [Nz_grid];
	T     = new double [Nz_grid];
	rho   = new double [Nz_grid];
	Pr    = new double [Nz_grid];
	c_eff = new double [Nz_grid]; // to be filled in getModalTrace(), depends on azimuth
	alpha = new double [Nz_grid];

	// Set up units of atmospheric profile object
	atm_profile->convert_altitude_units( Units::fromString( "m" ) );
	atm_profile->convert_property_units( "Z0", Units::fromString( "m" ) );
	atm_profile->convert_property_units( "U", Units::fromString( "m/s" ) );
	atm_profile->convert_property_units( "V", Units::fromString( "m/s" ) );
	atm_profile->convert_property_units( "T", Units::fromString( "K" ) );
	atm_profile->convert_property_units( "P", Units::fromString( "Pa" ) );
	atm_profile->convert_property_units( "RHO", Units::fromString( "kg/m3" ) );
  
	//
	// !!! ensure maxheight is less than the max height covered by the provided atm profile
	// may avoid some errors asociated with the code thinking that it goes above 
	// the max height when in fact the height may only differ from max height by
	// a rounding error. This should be revisited.
	// @todo revisit this
	// @todo add max_valid_height to AtmosphericProfile class
	
	if (maxheight >= atm_profile->get_maximum_altitude()) {
		maxheight = atm_profile->get_maximum_altitude() - 1e-6;
		cout << endl << "maxheight adjusted to: " << maxheight
			 << " m (max available in profile file)" << endl;
	}
  
	// fill and convert to SI units
	double dz       = (maxheight - z_min)/(Nz_grid - 1);	// the z-grid spacing
	//double z_min_km = z_min / 1000.0;
	//double dz_km    = dz / 1000.0;
  
	// Note: the rho, Pr, T, zw, mw are computed wrt ground level i.e.
	// the first value is at the ground level e.g. rho[0] = rho(z_min)
	// @todo make fill_vector( zvec ) methods in AtmosphericProfile()
	// @todo add underscores to internally calculated parameter keys
	atm_profile->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO", Units::fromString( "m/s" ) );
	atm_profile->calculate_wind_speed( "_WS_", "U", "V" );
	atm_profile->calculate_wind_direction( "_WD_", "U", "V" );
	atm_profile->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", freq );

	for (int i=0; i<Nz_grid; i++) {
		Hgt[i] = z_min + i*dz; // Hgt[0] = zground MSL
		rho[i] = atm_profile->get( "RHO", Hgt[i]);
		Pr[i]  = atm_profile->get( "P", Hgt[i] );
		T[i]   = atm_profile->get( "T", Hgt[i] );
		zw[i]  = atm_profile->get( "U", Hgt[i] );
		mw[i]  = atm_profile->get( "V", Hgt[i] );
		alpha[i]   = atm_profile->get( "_ALPHA_", Hgt[i] );
	}
}


// utility to print the parameters to the screen
void NCPA::ModeSolver::printParams() {
	printf(" Normal Modes Solver Parameters:\n");
	printf("                   freq : %g\n", freq);
	if (!Nby2Dprop) {
		printf("                azimuth : %g\n", azi);
	}
	else {
		printf("     azimuth_start (deg): %g\n", azi_min);
		printf("       azimuth_end (deg): %g\n", azi_max);
		printf("      azimuth_step (deg): %g\n", azi_step);
	}
	printf("                Nz_grid : %d\n", Nz_grid);
	printf("      z_min (meters MSL): %g\n", z_min);
	printf("      maxheight_km (MSL): %g\n", maxheight/1000.0);
	printf("   sourceheight_km (AGL): %g\n", sourceheight/1000.0);
	printf(" receiverheight_km (AGL): %g\n", receiverheight/1000.0);   
	printf("             Nrng_steps : %d\n", Nrng_steps);
	printf("            maxrange_km : %g\n", maxrange/1000.0); 
	printf("          gnd_imp_model : %s\n", gnd_imp_model.c_str());
	printf("Lamb wave boundary cond : %d\n", Lamb_wave_BC);
	printf("  SLEPc tolerance param : %g\n", tol);
	printf("    write_2D_TLoss flag : %d\n", write_2D_TLoss);
	printf("write_phase_speeds flag : %d\n", write_phase_speeds);
	printf("      write_speeds flag : %d\n", write_speeds);
	printf("  write_dispersion flag : %d\n", write_dispersion);
	printf("       write_modes flag : %d\n", write_modes);
	printf("         Nby2Dprop flag : %d\n", Nby2Dprop);
	printf("       turnoff_WKB flag : %d\n", turnoff_WKB);
	printf("    atmospheric profile : %s\n", atmosfile.c_str());
	if (!usrattfile.empty()) {
		printf("  User attenuation file : %s\n", usrattfile.c_str());
	}
	if (!modstartfile.empty()) {
		printf(" modal starter saved in : %s\n", modstartfile.c_str());
	}

	printf("       wvnum_filter_flg : %d\n", wvnum_filter_flg);
	if (wvnum_filter_flg==1) {
		printf("                  c_min : %g m/s\n", c_min);
		printf("                  c_max : %g m/s\n", c_max);
	}
}

*/







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
  
	FILE *tloss_1d    = fopen( filename_lossy.c_str(), "w");
	FILE *tloss_ll_1d = fopen( filename_lossless.c_str(), "w");
	// @todo check that files were opened properly

	for (i=0; i<n_r; i++) {
		r = (i+1)*dr;
		modal_sum_c    = 0.;
		modal_sum_c_ll = 0;
		modal_sum_i    = 0;      
		modal_sum_i_ll = 0;
      
		// Use the commented lines if the modes must be scaled by sqrt(rho)
		// @todo under what conditions will this be the case?  If never, we should remove entirely
		/*
		if (1) {
		//cout << sqrtrho_ratio << endl;
		for (m=0; m<select_modes; m++) {
		modal_sum_c    = modal_sum_c + v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio*exp(I*k_pert[m]*r)/sqrt(k_pert[m]);
		modal_sum_c_ll = modal_sum_c_ll + v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio*exp(I*real(k_pert[m])*r)/sqrt(real(k_pert[m]));
          
		modal_sum_i    = modal_sum_i + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio,2)*exp(-2*imag(k_pert[m])*r)/abs(k_pert[m]);
		modal_sum_i_ll = modal_sum_i_ll + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m]*sqrtrho_ratio,2)/real(k_pert[m]);
		}
		}
		*/

		// here we save the modal sum; the reduced pressure is: 
		//     p_red(r,z) =modal_sum/(4*pi*sqrt(rho(zs)))= p(r,z)/sqrt(rho(z))
		for (m=0; m<select_modes; m++) {    
			modal_sum_c    = modal_sum_c + v_s[n_zsrc][m] * v_s[n_zrcv][m]
					 * exp(I*k_pert[m]*r) / sqrt(k_pert[m]);
			modal_sum_c_ll = modal_sum_c_ll + v_s[n_zsrc][m] * v_s[n_zrcv][m]
					 * exp(I*real(k_pert[m])*r) / sqrt(real(k_pert[m]));
			modal_sum_i    = modal_sum_i + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2) 
					 * exp(-2*imag(k_pert[m])*r) / abs(k_pert[m]);
			modal_sum_i_ll = modal_sum_i_ll + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)
					 / real(k_pert[m]);
		}
      
		// no sqrt(rho[n_zrcv]/rho[n_zsrc]) factor
		if (1) {
			modal_sum_c    = expov8pi*modal_sum_c/sqrt(r);
			modal_sum_c_ll = expov8pi*modal_sum_c_ll/sqrt(r);
			modal_sum_i    = 4*PI*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
			modal_sum_i_ll = 4*PI*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		}
      
		// sqrt(rho[n_zrcv]/rho[n_zsrc]) factor added
		if (0) {
			modal_sum_c    = sqrtrho_ratio*expov8pi*modal_sum_c/sqrt(r);
			modal_sum_c_ll = sqrtrho_ratio*expov8pi*modal_sum_c_ll/sqrt(r);
			modal_sum_i    = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
			modal_sum_i_ll = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		}      

		fprintf(tloss_1d,"%f %20.12e %20.12e %20.12e\n", r/1000.0, real(modal_sum_c), 
			imag(modal_sum_c), modal_sum_i);
		fprintf(tloss_ll_1d,"%f %20.12e %20.12e %20.12e\n", r/1000.0, real(modal_sum_c_ll), 
			imag(modal_sum_c_ll), modal_sum_i_ll);
	}
	fclose(tloss_1d);
	fclose(tloss_ll_1d);
	printf("           file %s created\n", filename_lossy.c_str() );
	printf("           file %s created\n", filename_lossless.c_str() );
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
		tloss_ll_1d = fopen(filename_lossless.c_str(),"w");
	}
	else {  // append
		tloss_1d    = fopen(filename_lossy.c_str(),"a");
		tloss_ll_1d = fopen(filename_lossless.c_str(),"a");
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
			modal_sum_c_ll = modal_sum_c_ll 
				+ v_s[n_zsrc][m]*v_s[n_zrcv][m]*exp(I*real(k_pert[m])*r)/sqrt(real(k_pert[m]));
			modal_sum_i    = modal_sum_i    
				+ pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)*exp(-2*imag(k_pert[m])*r)/abs(k_pert[m]);
			modal_sum_i_ll = modal_sum_i_ll + pow(v_s[n_zsrc][m]*v_s[n_zrcv][m],2)/real(k_pert[m]);
		}
      
      
		// no sqrt(rho[n_zrcv]/rho[n_zsrc]) factor
		// @todo Need programmatic logic to determine which branch to take, or remove one
		if (1) {
			modal_sum_c    = expov8pi*modal_sum_c/sqrt(r);
			modal_sum_c_ll = expov8pi*modal_sum_c_ll/sqrt(r);
			modal_sum_i    = 4*PI*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
			modal_sum_i_ll = 4*PI*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		}
      
		// sqrt(rho[n_zrcv]/rho[n_zsrc]) factor added
		if (0) {
			modal_sum_c    = sqrtrho_ratio*expov8pi*modal_sum_c/sqrt(r);
			modal_sum_c_ll = sqrtrho_ratio*expov8pi*modal_sum_c_ll/sqrt(r);
			modal_sum_i    = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i)*sqrt(1./8./PI/r);
			modal_sum_i_ll = 4*PI*sqrtrho_ratio*sqrt(modal_sum_i_ll)*sqrt(1./8./PI/r);
		}       
      
		fprintf(tloss_1d,"%10.3f %8.3f %20.12e %20.12e %20.12e\n", r/1000.0, 
			azimuth, real(modal_sum_c), imag(modal_sum_c), modal_sum_i);
		fprintf(tloss_ll_1d,"%10.3f %8.3f %20.12e %20.12e %20.12e\n", r/1000.0, 
			azimuth, real(modal_sum_c_ll), imag(modal_sum_c_ll), modal_sum_i_ll);
	}
	fprintf(tloss_1d, "\n");
	fprintf(tloss_ll_1d, "\n");
  
	fclose(tloss_1d);
	fclose(tloss_ll_1d);
	return 0;
}

// 20150603 DV added code to account for sqrt(rho_rcv/rho_src) if so chosen
int NCPA::ModeSolver::getTLoss2D(int nz, int select_modes, double dz, int n_r, double dr, 
	double z_src, double *rho, complex<double> *k_pert, double **v_s, 
	std::string filename_lossy) {

	// !!! note the modes assumed here are the modes in Oc. Acoust. divided by sqrt(rho(z))
	// Thus the formula for (physical) pressure is:
	// p(r,z) = i*exp(-i*pi/4)/sqrt(8*pi*r)*1/rho(zs)* sum[ v_s(zs)*sqrt(rho(zs))*v_s(z)*sqrt(rho(z))
	// *exp(ikr)/sqrt(k) ]
	// We usually save the TL corresponding to the reduced pressure: p_red(r,z) = p(r,z)/sqrt(rho(z))
  
	int i, j, m, stepj;
	int n_zsrc = (int) ceil(z_src/dz);
	double r, z, sqrtrhoj, rho_atzsrc;
	complex<double> modal_sum;
	complex<double> I (0.0, 1.0);
	
	// the 4*PI factor ensures that the modal sum below ends up being the actual TL
	complex<double> expov8pi = 4*PI*I*exp(-I*PI*0.25)/sqrt(8.0*PI); 
  
	rho_atzsrc = rho[n_zsrc];

	stepj = nz/500; // controls the vertical sampling of 2D data saved 
	if (stepj==0) {
		stepj = 10;	// ensure it's never 0; necessary for the loop below
	}

	FILE *tloss_2d = fopen(filename_lossy.c_str(),"w");

	for (i=0; i<n_r; i++) {
		r = (i+1)*dr;
		for (j=0; j<nz; j=j+stepj) {
			z = (j)*dz;
			sqrtrhoj = sqrt(rho[j]);
			modal_sum = 0.;
          
			// Use the commented lines if the modes must be scaled by sqrt(rho)
			// @todo do we need these?
			/*
			if (1) {
			for (m=0; m<select_modes; m++) {
			modal_sum = modal_sum + v_s[n_zsrc][m]*v_s[j][m]*sqrtrhoj*exp(I*k_pert[m]*r)/sqrt(k_pert[m]);
			}
			modal_sum = expov8pi*modal_sum/sqrt(r*rho_atzsrc);
			}
			*/
                  
			for (m=0; m<select_modes; m++) {
				modal_sum = modal_sum + v_s[n_zsrc][m]*v_s[j][m]*exp(I*k_pert[m]*r)/sqrt(k_pert[m]);
			}
			modal_sum = expov8pi*modal_sum/sqrt(r); // no sqrt(rho[n_zrcv]/rho[n_zsrc]) factor

			// @todo do we need this?
			if (0) {
				// sqrt(rho[n_zrcv]/rho[n_zsrc]) factor added
				modal_sum = sqrtrhoj/sqrt(rho_atzsrc)*modal_sum/sqrt(r);
			}           

			fprintf(tloss_2d,"%f %f %15.8e %15.8e\n", r/1000.0, z/1000.0, 
				real(modal_sum), imag(modal_sum));
		}
		fprintf(tloss_2d,"\n");
	}
	fclose(tloss_2d);
	printf("           file %s created\n", filename_lossy.c_str() );
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
 

int NCPA::ModeSolver::writePhaseSpeeds(int select_modes, double freq, complex<double> *k_pert)
{
	int j;
	FILE *phasespeeds= fopen(tag_filename("phasespeeds.nm").c_str(), "w");
	// @todo check to see if file is opened
	for (j=0; j< select_modes; j++) {
		fprintf(phasespeeds, "%d %f %15.8e\n", j, (2*PI*freq)/real(k_pert[j]), imag(k_pert[j]));
	}
	fclose(phasespeeds);
	printf("           file %s created\n",tag_filename("phasespeeds.nm").c_str());
	return 0;
}



int NCPA::ModeSolver::writeEigenFunctions(int nz, int select_modes, double dz, double **v_s)
{
	char mode_output[40];
	int j, n;
	double chk;
	double dz_km = dz/1000.0;

	for (j=0; j<select_modes; j++) {
		sprintf(mode_output,"mode_%d.nm", j);
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
	printf("           files mode_<mode_number> created (%d in total)\n", select_modes);  
	return 0;
}

std::string NCPA::ModeSolver::tag_filename( std::string basename ) {
	return user_tag + basename;
}