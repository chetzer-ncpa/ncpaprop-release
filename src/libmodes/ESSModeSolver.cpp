#include "slepceps.h"
#include "slepcst.h"
#include <complex>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "modes.h"
#include "EigenEngine.h"
#include "Atmosphere1D.h"
#include "util.h"

#define MAX_MODES 4000 

using namespace NCPA;
using namespace std;


NCPA::ESSModeSolver::ESSModeSolver( NCPA::ParameterSet *param, NCPA::Atmosphere1D *atm_profile )
{
	setParams( param, atm_profile );           
}

NCPA::ESSModeSolver::~ESSModeSolver() {}

void NCPA::ESSModeSolver::setParams( NCPA::ParameterSet *param, NCPA::Atmosphere1D *atm_prof )
{		

	// obtain the parameter values from the user's options
	// @todo add units to input scalar quantities
	atmosfile 			= param->getString( "atmosfile" );
	gnd_imp_model 		= param->getString( "ground_impedence_model" );
	usrattfile 			= param->getString( "use_attn_file" );
	modstartfile 		= param->getString( "modal_starter_file" );
  	z_min 				= param->getFloat( "zground_km" ) * 1000.0;    // meters
  	maxrange 			= param->getFloat( "maxrange_km" ) * 1000.0;
  	maxheight 			= param->getFloat( "maxheight_km" ) * 1000.0;      // @todo fix elsewhere that m is required
  	sourceheight 		= param->getFloat( "sourceheight_km" ) * 1000.0;
  	receiverheight 		= param->getFloat( "receiverheight_km" ) * 1000.0;
  	tol 				= 1.0e-8;
  	Nz_grid 			= param->getInteger( "Nz_grid" );
  	Nrng_steps 			= param->getInteger( "Nrng_steps" );
  	Lamb_wave_BC 		= param->getBool( "Lamb_wave_BC" );
  	write_2D_TLoss  	= param->getBool( "write_2d_tloss" );
  	write_atmosphere    = param->getBool( "write_atm_profile" );
  	write_phase_speeds 	= param->getBool( "write_phase_speeds" );
  	write_speeds 		= param->getBool( "write_speeds" );
  	write_modes 		= param->getBool( "write_modes" );
  	//write_dispersion 	= param->getBool( "write_dispersion" );
  	dispersion_file     = param->getString( "dispersion_file" );
  	append_dispersion_file
  						= param->getBool( "append_dispersion_file" );
  	Nby2Dprop 			= param->getBool( "multiprop" );
  	turnoff_WKB 		= param->getBool( "turnoff_WKB" );
  	z_min_specified     = param->wasFound( "zground_km" );
  	user_tag			= param->getString( "filetag" );
  	if (user_tag.size() > 0) {
  		user_tag += ".";
  	}

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
			       || write_modes || (dispersion_file.size() > 0)) {
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

	// frequencies
	if (param->getBool( "broadband" ) ) {
		double f_min, f_step, f_max;
      	f_min = param->getFloat( "f_min" );
      	f_step = param->getFloat( "f_step" );
      	f_max = param->getFloat( "f_max" );

      	// sanity checks
      	if (f_min >= f_max) {
      		throw std::runtime_error( "f_min must be less than f_max!" );
      	}

      	Nfreq = (double)(std::floor((f_max - f_min) / f_step)) + 1;
      	f_vec = new double[ Nfreq ];
      	std::memset( f_vec, 0, Nfreq * sizeof( double ) );
      	for (int fi = 0; fi < Nfreq; fi++) {
      		f_vec[ fi ] = f_min + ((double)fi) * f_step;
      	}
	} else {
		freq = param->getFloat( "freq" );
		f_vec = new double[ 1 ];
		f_vec[ 0 ] = freq;
		Nfreq = 1;
	}

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
	atm_profile->convert_property_units( "U", Units::fromString( "m/s" ) );
	atm_profile->convert_property_units( "V", Units::fromString( "m/s" ) );
	atm_profile->convert_property_units( "T", Units::fromString( "K" ) );
	atm_profile->convert_property_units( "P", Units::fromString( "Pa" ) );
	atm_profile->convert_property_units( "RHO", Units::fromString( "kg/m3" ) );
	if (atm_profile->contains_scalar("Z0")) {
		atm_profile->convert_property_units( "Z0", Units::fromString( "m" ) );
		double profile_z0 = atm_profile->get( "Z0" );
		if (z_min < profile_z0) {
			std::cout << "Adjusting minimum altitude to profile ground height of "
			          << profile_z0 << " m" << std::endl;
			z_min = profile_z0;
		}
	}
  
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
	if (atm_profile->contains_vector( "C0" ) ) {
		std::cout << "Using user-supplied static sound speed" << std::endl;
		atm_profile->copy_vector_property( "C0", "_C0_" );
		atm_profile->convert_property_units( "_C0_", Units::fromString( "m/s" ) );
	} else {
		if (atm_profile->contains_vector( "P" ) && atm_profile->contains_vector("RHO")) {
			std::cout << "Calculating sound speed from pressure and density." << std::endl;
			atm_profile->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO", 
				Units::fromString( "m/s" ) );
		} else {
			std::cout << "Calculating sound speed from temperature" << std::endl;
			atm_profile->calculate_sound_speed_from_temperature( "_C0_", "T", 
				Units::fromString( "m/s" ) );
		}
	}
	atm_profile->calculate_wind_speed( "_WS_", "U", "V" );
	atm_profile->calculate_wind_direction( "_WD_", "U", "V" );
	if (usrattfile.empty()) {
		atm_profile->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", freq );
	} else {
		atm_profile->read_attenuation_from_file( "_ALPHA_", usrattfile );
	}
	

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
void NCPA::ESSModeSolver::printParams() {
	printf(" Effective Sound Speed Normal Modes Solver Parameters:\n");
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
	printf("    write_2d_tloss flag : %d\n", write_2D_TLoss);
	printf("write_phase_speeds flag : %d\n", write_phase_speeds);
	printf("      write_speeds flag : %d\n", write_speeds);
	printf("       write_modes flag : %d\n", write_modes);
	printf("         multiprop flag : %d\n", Nby2Dprop);
	printf("       turnoff_WKB flag : %d\n", turnoff_WKB);
	printf("    atmospheric profile : %s\n", atmosfile.c_str());
	if (!dispersion_file.empty()) {
		printf("        dispersion file : %s\n", dispersion_file.c_str());
		printf(" Dispersion file append : %d\n", append_dispersion_file );
	}
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

// int NCPA::ESSModeSolver::doESSSLEPcCalculation( double *diag, double dz, double *k_min, 
// 	double *k_max, PetscInt *nconv, double *k2, double **v ) {

// 	//
// 	// Declarations related to Slepc computations
// 	//
// 	Mat            A;           // problem matrix
// 	EPS            eps;         // eigenproblem solver context
// 	ST             stx;
// 	KSP            kspx;
// 	PC             pcx;
// 	EPSType        type;        // CHH 191022: removed const qualifier
// 	PetscReal      re, im;
// 	PetscScalar    kr, ki, *xr_;
// 	Vec            xr, xi;
// 	PetscInt       Istart, Iend, col[3], its, maxit;
// 	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
// 	PetscScalar    value[3];	
// 	PetscErrorCode ierr;
// 	PetscMPIInt    rank, size;
// 	int i, j, nev;
// 	double h2 = dz * dz;


// 	// Initialize Slepc
// 	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
// 	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);  

// 	// Create the matrix A to use in the eigensystem problem: Ak=kx
// 	ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
// 	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nz_grid,Nz_grid); CHKERRQ(ierr);
// 	ierr = MatSetFromOptions(A); CHKERRQ(ierr);

// 	// the following Preallocation call is needed in PETSc version 3.3
// 	ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL); CHKERRQ(ierr);
// 	// or use: ierr = MatSetUp(A); 

// 	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 	Compute the operator matrix that defines the eigensystem, Ax=kx
// 	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// 	// Make matrix A 
// 	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
// 	if (Istart==0) 
// 		FirstBlock=PETSC_TRUE;
// 	if (Iend==Nz_grid) 		/* @todo check if should be Nz_grid-1 */
// 		LastBlock=PETSC_TRUE;    

// 	value[0]=1.0/h2; 
// 	value[2]=1.0/h2;
// 	for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
// 		value[1] = -2.0/h2 + diag[i];
// 		col[0]=i-1; 
// 		col[1]=i; 
// 		col[2]=i+1;
// 		ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES); CHKERRQ(ierr);
// 	}
// 	if (LastBlock) {
// 		i=Nz_grid-1; 
// 		col[0]=Nz_grid-2; 
// 		col[1]=Nz_grid-1;
// 		ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
// 	}
// 	if (FirstBlock) {
// 		i=0; 
// 		col[0]=0; 
// 		col[1]=1; 
// 		value[0]=-2.0/h2 + diag[0]; 
// 		value[1]=1.0/h2;
// 		ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
// 	}

// 	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
// 	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

// 	// CHH 191022: MatGetVecs() is deprecated, changed to MatCreateVecs()
// 	ierr = MatCreateVecs(A,PETSC_NULL,&xr); CHKERRQ(ierr);
// 	ierr = MatCreateVecs(A,PETSC_NULL,&xi); CHKERRQ(ierr);

// 	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 	Create the eigensolver and set various options
// 	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 	/* 
// 	Create eigensolver context
// 	*/
// 	ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

// 	/* 
// 	Set operators. In this case, it is a standard eigenvalue problem
// 	*/
// 	ierr = EPSSetOperators(eps,A,PETSC_NULL); CHKERRQ(ierr);
// 	ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRQ(ierr);

// 	/*
// 	Set solver parameters at runtime
// 	*/
// 	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
// 	ierr = EPSSetType(eps,"krylovschur"); CHKERRQ(ierr);
// 	ierr = EPSSetDimensions(eps,10,PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr); // leaving this line in speeds up the code; better if this is computed in chunks of 10? - consult Slepc manual
// 	ierr = EPSSetTolerances(eps,tol,PETSC_DECIDE); CHKERRQ(ierr);

// 	ierr = EPSGetST(eps,&stx); CHKERRQ(ierr);
// 	ierr = STGetKSP(stx,&kspx); CHKERRQ(ierr);
// 	ierr = KSPGetPC(kspx,&pcx); CHKERRQ(ierr);
// 	ierr = STSetType(stx,"sinvert"); CHKERRQ(ierr);
// 	ierr = KSPSetType(kspx,"preonly");
// 	ierr = PCSetType(pcx,"cholesky");
// 	ierr = EPSSetInterval(eps,pow(*k_min,2),pow(*k_max,2)); CHKERRQ(ierr);
// 	ierr = EPSSetWhichEigenpairs(eps,EPS_ALL); CHKERRQ(ierr);
	
// 	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 	Solve the eigensystem
// 	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// 	ierr = EPSSolve(eps);CHKERRQ(ierr);
// 	/*
// 	Optional: Get some information from the solver and display it
// 	*/
// 	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
// 	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
// 	ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
// 	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	
// 	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// 	Display solution and clean up
// 	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
// 	/* 
// 	Get number of converged approximate eigenpairs
// 	*/
// 	ierr = EPSGetConverged(eps,nconv);CHKERRQ(ierr);
	
// 	if ((*nconv)>0) {
// 		for (i=0;i<(*nconv);i++) {
// 			ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
// #if defined(PETSC_USE_COMPLEX)
// 			re = PetscRealPart(kr);
// 			im = PetscImaginaryPart(kr);
// #else
// 			re = kr;
// 			im = ki;
// #endif 
// 			k2[(*nconv)-i-1] = re; // proper count of modes
// 			ierr = VecGetArray(xr,&xr_);CHKERRQ(ierr);
// 			for (j = 0; j < Nz_grid; j++) {
// 				v[j][(*nconv)-i-1] = xr_[j]/sqrt(dz); //per Slepc the 2-norm of xr_ is=1; we need sum(v^2)*dz=1 hence the scaling xr_/sqrt(dz)
// 			}
// 			ierr = VecRestoreArray(xr,&xr_);CHKERRQ(ierr);
// 		}
// 	}

// 	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
// 	ierr = MatDestroy(&A);  CHKERRQ(ierr);
// 	ierr = VecDestroy(&xr); CHKERRQ(ierr);
// 	ierr = VecDestroy(&xi); CHKERRQ(ierr); 

	

// 	return 0;
// }



int NCPA::ESSModeSolver::solve() {
	//
	// Declarations related to Slepc computations
	//
	/*
	Mat            A;           // problem matrix
	EPS            eps;         // eigenproblem solver context
	ST             stx;
	KSP            kspx;
	PC             pcx;
	EPSType        type;        // CHH 191022: removed const qualifier
	PetscReal      re, im;
	PetscScalar    kr, ki, *xr_;
	Vec            xr, xi;
	PetscInt       Istart, Iend, col[3], its, maxit, nconv;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar    value[3];	
	PetscErrorCode ierr;
	PetscMPIInt    rank, size;
	*/
	PetscInt nconv;
	PetscErrorCode ierr;

	int    i, select_modes, nev, it, fi;
	double dz, admittance, rng_step, z_min_km;  // , h2;
	double k_min, k_max;			
	double *diag, *k2, *k_s, **v, **v_s;
	complex<double> *k_pert;

	diag   = new double [Nz_grid];
  
	k2     = new double [MAX_MODES];
	k_s    = new double [MAX_MODES];
	k_pert = new complex<double> [MAX_MODES];
	v      = dmatrix(Nz_grid,MAX_MODES);
	v_s    = dmatrix(Nz_grid,MAX_MODES); 			

	//nev   = 0;
	k_min = 0; 
	k_max = 0;  	

	rng_step = maxrange/Nrng_steps;         		// range step [meters]
	dz       = (maxheight - z_min)/(Nz_grid - 1);	// the z-grid spacing
	//h2       = dz*dz;
	z_min_km = z_min / 1000.0;

	SlepcInitialize(PETSC_NULL,PETSC_NULL,(char*)0,PETSC_NULL); /* @todo move out of loop? */

	for (fi = 0; fi < Nfreq; fi++) {
		freq = f_vec[ fi ];

		//
		// loop over azimuths (if not (N by 2D) it's only one azimuth)
		//
		for (it=0; it<Naz; it++) {
	  
			azi = azi_min + it*azi_step; // degrees (not radians)
			cout << endl << "Now processing azimuth = " << azi << " (" << it+1 << " of " 
				 << Naz << ")" << endl;

			atm_profile->calculate_wind_component("_WC_", "_WS_", "_WD_", azi );
			atm_profile->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );
			atm_profile->get_property_vector( "_CEFF_", c_eff );
			atm_profile->add_property( "_AZ_", azi, NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );

			//
			// ground impedance model
			//
			// at the ground the BC is: Phi' = (a - 1/2*dln(rho)/dz)*Phi
			// for a rigid ground a=0; and the BC is the Lamb Wave BC:
			// admittance = -1/2*dln(rho)/dz
			//  
			admittance = 0.0;   // default
			if ( gnd_imp_model.compare( "rigid" ) == 0) {
				// Rigid ground, check for Lamb wave BC
				if (Lamb_wave_BC) {
					admittance = -atm_profile->get_first_derivative( "RHO", z_min ) 
						/ atm_profile->get( "RHO", z_min ) / 2.0;
					cout << "Admittance = " << admittance << endl;
				} else {
					admittance = 0.0;
				}
			} else {
				throw invalid_argument( 
					"This ground impedance model is not implemented yet: " + gnd_imp_model );
			}

			//
			// Get the main diagonal and the number of modes
			//		
			i = getModalTrace(Nz_grid, z_min, sourceheight, receiverheight, dz, atm_profile, 
					admittance, freq, azi, diag, &k_min, &k_max, turnoff_WKB, c_eff);

			// if wavenumber filtering is on, redefine k_min, k_max
			if (wvnum_filter_flg) {
				k_min = 2 * PI * freq / c_max;
				k_max = 2 * PI * freq / c_min;
			}

			i = getNumberOfModes(Nz_grid,dz,diag,k_min,k_max,&nev);

			printf ("______________________________________________________________________\n\n");
			printf (" -> Normal mode solution at %5.3f Hz and %5.2f deg (%d modes)...\n", freq, azi, nev);
			printf (" -> Discrete spectrum: %5.2f m/s to %5.2f m/s\n", 2*PI*freq/k_max, 2*PI*freq/k_min);
	    
	    	i = NCPA::EigenEngine::doESSCalculation( diag, Nz_grid, dz, tol,
	    		&k_min, &k_max, &nconv, k2, v );

			// select modes and do perturbation
			doSelect(Nz_grid,nconv,k_min,k_max,k2,v,k_s,v_s,&select_modes);
			doPerturb(Nz_grid, z_min, dz, select_modes, freq, k_s, v_s, alpha, k_pert);
			
			//
			// Output data  
			//
			if (Nby2Dprop) { // if (N by 2D is requested)
				getTLoss1DNx2(azi, select_modes, dz, Nrng_steps, rng_step, sourceheight, 
					      receiverheight, rho, k_pert, v_s, Nby2Dprop, it,
					      tag_filename("Nby2D_tloss_1d.nm"),
					      tag_filename("Nby2D_tloss_1d.lossless.nm") );
			}
			else {
				cout << "Writing to file: 1D transmission loss at the ground..." << endl;
				getTLoss1D(select_modes, dz, Nrng_steps, rng_step, sourceheight, receiverheight, rho, 
					   k_pert, v_s, tag_filename("tloss_1d.nm"),
					   tag_filename("tloss_1d.lossless.nm") );

				if ( !modstartfile.empty() ) {
					cout << "Writing to file: modal starter" << endl;
					
					// Modal starter - DV 20151014
					// Modification to apply the sqrt(k0) factor to agree with Jelle's getModalStarter; 
					// this in turn will make 'pape' agree with modess output
					getModalStarter(Nz_grid, select_modes, dz, freq, sourceheight, receiverheight, 
							rho, k_pert, v_s, tag_filename(modstartfile) );
				}

				if (write_2D_TLoss) {
					cout << "Writing to file: 2D transmission loss...\n";
					getTLoss2D(Nz_grid,select_modes,dz,Nrng_steps,rng_step,sourceheight,rho, 
						   k_pert,v_s, tag_filename("tloss_2d.nm") );
				}

				if (write_phase_speeds) {
					cout << "Writing to file: phase speeds...\n";
					writePhaseSpeeds(select_modes,freq,k_pert);
				}

				if (write_modes) {
					cout << "Writing to file: the modes and the phase and group speeds...\n";
					writeEigenFunctions(Nz_grid,select_modes,dz,v_s);
					writePhaseAndGroupSpeeds(Nz_grid, dz, select_modes, freq, k_pert, v_s, c_eff);
				}
	        
				if ((write_speeds) &&  !(write_modes)) {
					cout << "Writing to file: the modal phase speeds and the group speeds...\n";
					writePhaseAndGroupSpeeds(Nz_grid, dz, select_modes, freq, k_pert, v_s, c_eff);
				}        

				if (!dispersion_file.empty()) {
					printf ("Writing to file %s: dispersion at freq = %8.3f Hz...\n",
						dispersion_file.c_str(), freq);
					FILE *dispersionfile;
					if (append_dispersion_file) {
						dispersionfile = fopen(dispersion_file.c_str(),"a");
					} else {
						dispersionfile = fopen(dispersion_file.c_str(),"w");
					}
					writeDispersion(dispersionfile,select_modes,dz,sourceheight,receiverheight,
						freq,k_pert,v_s, rho);
					fclose(dispersionfile);
				}

				if (write_atmosphere) {
					std::cout << "Writing atmosphere to "
						<< tag_filename("atm_profile.nm") << std::endl;
					std::vector<std::string> keylist;
					keylist.push_back( "U" );
					keylist.push_back( "V" );
					keylist.push_back( "T" );
					keylist.push_back( "RHO" );
					keylist.push_back( "P" );
					keylist.push_back( "_C0_" );
					keylist.push_back( "_CEFF_" );
					std::ofstream atmout( tag_filename( "atm_profile.nm" ) );
					atm_profile->print_atmosphere( keylist, "Z", atmout );
					atmout.close();
				}
			}
	    
			// Clean up azimuth-specific atmospheric properties before starting next run
			atm_profile->remove_property( "_WC_" );
			atm_profile->remove_property( "_CEFF_" );
			atm_profile->remove_property( "_AZ_" );
			
	  
		} // end loop by azimuths

	} // end loop by frequencies
  
	// Finalize Slepc
	ierr = SlepcFinalize();CHKERRQ(ierr);

  
	// free the rest of locally dynamically allocated space
	delete[] diag;
	delete[] k2;
	delete[] k_s;
	delete[] k_pert;
	free_dmatrix(v, Nz_grid, MAX_MODES);
	free_dmatrix(v_s, Nz_grid, MAX_MODES);
  
	return 0;
} // end of computeModes()






// DV 20130827: new getModalTrace() to output the effective sound speed 
int NCPA::ESSModeSolver::getModalTrace( int nz, double z_min, double sourceheight, double receiverheight, 
	double dz, Atmosphere1D *p, double admittance, double freq, double azi, double *diag,
	double *k_min, double *k_max, bool turnoff_WKB, double *ceffz) {
	
	// use atmospherics input for the trace of the matrix.
	// the vector diag can be used to solve the modal problem
	// also returns the bounds on the wavenumber spectrum, [k_min,k_max]
	int    i, top;
	double azi_rad, z_km, z_min_km, dz_km, omega, gamma, bnd_cnd; 
	double cz, windz, ceffmin, ceffmax, ceff_grnd, cefftop; 
	double kk, dkk, k_eff, k_gnd, k_max_full, wkbIntegral, wkbTerm;
	z_min_km = z_min / 1000.0;
	dz_km    = dz / 1000.0;
	omega    = 2*PI*freq;
  
	azi_rad  = NCPA::Units::convert( p->get( "_AZ_" ), NCPA::UNITS_ANGLE_DEGREES, UNITS_ANGLE_RADIANS );
  
	gamma = 1.4;  
	
	z_km      = z_min_km;
	
	cz        = p->get( "_C0_", z_min );
	windz     = p->get( "_WC_", z_min );
	ceff_grnd = p->get( "_CEFF_", z_min );
	ceffmin   = ceff_grnd;  // in m/s; initialize ceffmin
	ceffmax   = ceffmin;    // initialize ceffmax   
	for (i=0; i<nz; i++) {	     
		ceffz[ i ] = p->get( "_CEFF_", Hgt[ i ] );
		// we neglect rho_factor - it does not make sense to have in this approximation
		//rho_factor is 1/2*rho_0"/rho_0 - 3/4*(rho_0')^2/rho_0^2
		diag[i] = pow( omega / ceffz[ i ], 2 );

		if (ceffz[i] < ceffmin)
			ceffmin = ceffz[i];  // approximation to find minimum sound speed in problem		   		
		if (ceffz[i] > ceffmax)
			ceffmax = ceffz[i];

		z_km += dz_km;		
	}
  
	bnd_cnd = (1./(dz*admittance+1))/(pow(dz,2)); // bnd cnd assuming centered fd
	diag[0] = bnd_cnd + diag[0]; 
  
	// use WKB trick for ground to ground propagation.
	bool used_WKB = false;
	if ((fabs(sourceheight)<1.0e-3) && (fabs(receiverheight)<1.0e-3) && (!turnoff_WKB)) {
		//
		// WKB trick for ground to ground propagation. 
		// Cut off lower phasespeed (highest wavenumber) where tunneling 
		// is insignificant (freq. dependent)
		//
		used_WKB = true;
		k_max_full = omega/ceffmin; 
		k_gnd      = omega/ceff_grnd;
		dkk        = (pow(k_max_full,2) - pow(k_gnd,2)) / 100.0;
    
		//
		// dv: when ceffmin is at the ground dkk can be very small but non-zero (rounding errors?)
		// in that case we need to skip the next {for loop}; otherwise it will be executed
		// with the small increment dkk
		//
		kk = pow(k_gnd,2);   //initialization of kk in case the loop is skipped
		if (dkk >1.0e-10) {  // do this loop only if dkk is nonzero (significantly)
			for (kk = pow(k_gnd,2); kk < pow(k_max_full,2); kk=kk+dkk) {
				wkbIntegral = 0.0;
				wkbTerm     = 1.0;  
				i           = 0;
				z_km        = z_min_km;
				while (wkbTerm > dkk) {
					k_eff       = omega/ceffz[i];
					wkbTerm     = abs(kk - pow(k_eff,2));
					wkbIntegral = wkbIntegral + dz*sqrt(wkbTerm); // dz should be in meters					
					i++;
					z_km += dz_km;
				} 

				if (wkbIntegral >= 10.0) {
					printf("\nWKB fix: new phasevelocity minimum= %6.2f m/s (was %6.2f m/s)\n", \
						omega/sqrt(kk), omega/k_max_full); 
					break;
				}
			}
		}  

		*k_max = sqrt(kk);  // use this for ground-to-ground 1D Tloss 
		//(uses WKB trick to include only non-vanishing modes at the ground)			  
		// *k_max = k_max_full; 
	} else { // not ground-to-ground propagation
		*k_max = omega/ceffmin; // same as k_max_full
	}  

	top     = nz - ((int) nz/10);
	z_km    = z_min_km + (top+1)*dz_km;  
	cz      = p->get( "_C0_", Hgt[ top+1 ] );
	windz   = p->get( "_WC_", Hgt[ top+1 ] );
	cefftop = p->get( "_CEFF_", Hgt[ top+1 ] );
	*k_min  = omega/cefftop;

	if (used_WKB && (*k_max < *k_min)) {
		std::cout << "Calculated k_max less than k_min, turning off WKB"
				  << std::endl;
		turnoff_WKB = true;
		*k_max = omega/ceffmin;
	}

	// optional save ceff
	// @todo add flag to turn on/off
	if (0) {
		double *target, *zvec;
		size_t nz = p->nz();
		target = new double[ nz ];
		zvec   = new double[ nz ];
		p->get_property_vector( "_CEFF_", target );
		p->get_altitude_vector( zvec );
		
		FILE *fp = fopen("ceff.nm", "w");    
		for ( size_t ii = 0; ii < nz; ii++ ) {     
			z_km = zvec[ ii ];
			fprintf(fp, "%8.3f %15.6e\n", zvec[ ii ], target[ ii ]);
		}
		fclose(fp);
		printf("ceff saved in ceff.nm\n");
		delete [] target;
		delete [] zvec;
	}

	return 0;
}

int NCPA::ESSModeSolver::doSelect(int nz, int n_modes, double k_min, double k_max, 
	double *k2, double **v, double *k_s, double **v_s, int *select_modes)
{
	int cnt = -1;
	int i, j;
	double k;
	for (j=0; j<n_modes; j++) {
		k = sqrt(k2[j]);
		if ((k >= k_min) && (k <= k_max)) {
			cnt++;
			for (i=0; i<nz; i++) {
				v_s[i][cnt] = v[i][j];
			}
			k_s[cnt] = k;
		}
	}
	*select_modes = cnt + 1; //+1 because cnt started from zero
	return 0;
}

// Modal starter - DV 20151014
// Modification to apply the sqrt(k0) factor to agree with Jelle's getModalStarter; 
// this in turn will make 'pape' agree with modess output
void NCPA::ESSModeSolver::getModalStarter(int nz, int select_modes, double dz, double freq,  
	double z_src, double z_rcv, double *rho, complex<double> *k_pert, double **v_s, 
	string modstartfile) {

	// computes the modal starter field useful to ingest in a subsequent PE run.
	// The formula for the modal starter as in Ocean Acoustics, 1994 ed. eq. 6.72, pg. 361:
	// where Psi_m(z) are the modes defined in (Ocean Acoustics, 1994 ed. pg. 274.)
	// Eq. 6.72 is the "normalized" modal field which in this case means
	// that it is divided by (1/(4*PI)) i.e. multiplied by (4*PI)
	// the abs(pressure of point source in free space).
	// We do not need this normalization so we'll use formula 6.72 / (4*PI)
  
  
	// Psi_starter(z) = sqrt(2*PI)/(4*PI)*sqrt(rho(z))/sqrt(rho(zs))*sum( Vm(z)*Vm(zs)/sqrt(k_m) )
	// where Vm are the modes computed in this code. Note that Vm(z) = Psi_m(z)/sqrt(rho(z)), 
	// 

	int j, m;
	int n_zsrc = (int) ceil(z_src/dz);
	double z;
	complex<double> modal_sum;
  
	double k0 = (2*PI*freq/340.0); // reference wavenumber
	int z_cnd = int ((340.0/freq)/10/dz);
	FILE *mstfile = fopen(modstartfile.c_str(),"w");

	for (j=0; j<nz; j=j+z_cnd) {
		z = (j)*dz;
		modal_sum = 0.0;
		
		// Use the commented lines if the modes must be scaled by sqrt(rho)
		for (m=0; m<select_modes; m++) {
			modal_sum = modal_sum + v_s[n_zsrc][m]*v_s[j][m]/sqrt(real(k_pert[m]));
		}     
      
		//modal_sum = factor/twoPI*sqrtTwo*sqrtrhoj*modal_sum*sqrt(k0); // modal starter at particular z
      
		modal_sum = modal_sum * PI*sqrt(k0); // Jelle - needs the sqrt(k0)
		//modal_sum = factor/twoPI*sqrtTwo*sqrtrhoj*modal_sum; // modal starter at particular z
               
		//for (m=0; m<select_modes; m++) {
		//    modal_sum = modal_sum + v_s[n_zsrc][m]*v_s[j][m]/sqrt(k_pert[m]);
		//}
		//modal_sum = sqrt2PI/rhos*modal_sum; // modal starter at particular z
		//modal_sum = factor/rhos*modal_sum; // modal starter at particular z
      
      
		fprintf(mstfile,"%10.3f   %16.12e   %16.12e\n", z/1000.0, real(modal_sum), imag(modal_sum));
	}
	fprintf(mstfile,"\n");
	fclose(mstfile);
	printf("           file %s created\n", modstartfile.c_str());  

}

// compute/save the group and phase speeds
int NCPA::ESSModeSolver::writePhaseAndGroupSpeeds(int nz, double dz, int select_modes, double freq, 
	complex<double> *k_pert, double **v_s, double *c_eff) {
	int j, n;
	double v_phase, v_group;
	double omega = 2*PI*freq;
  
	FILE *speeds= fopen(tag_filename("speeds.nm").c_str(), "w");
	for (j=0; j< select_modes; j++) {
		v_phase = omega/real(k_pert[j]);
      
		// compute the group speed: vg = v_phase*int_0^z_max Psi_j^2/c_eff^2 dz_km
		v_group = 0.0;
		for (n=0; n<nz; n++) {
			v_group = v_group + v_s[n][j]*v_s[n][j]/(c_eff[n]*c_eff[n]);   
		}
		v_group = v_group*v_phase*dz;
		v_group = 1.0/v_group;
      
		fprintf(speeds, "%4d %9.3f %9.3f %15.8e\n", j+1, v_phase, v_group, imag(k_pert[j]));    
	}
	fclose(speeds);
	printf("           Phase and group speeds saved in file %s.\n",
		tag_filename("speeds.nm").c_str());
	return 0;
}
