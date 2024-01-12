#include "slepceps.h"
#include "slepcst.h"
#include <complex>
#include <stdexcept>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "modes.h"
#include "EigenEngine.h"
#include "Atmosphere1D.h"
#include "util.h"

#define MAX_MODES 4000 

using namespace NCPA;
using namespace std;


NCPA::WModeSolver::WModeSolver( NCPA::ParameterSet *param, NCPA::Atmosphere1D *atm_profile )
{
	setParams( param, atm_profile );           
}

NCPA::WModeSolver::~WModeSolver() {}

void NCPA::WModeSolver::setParams( NCPA::ParameterSet *param, NCPA::Atmosphere1D *atm_prof )
{		

	// obtain the parameter values from the user's options
  // @todo add units to input scalar quantities
	atmosfile 			= param->getString( "atmosfile" );
	gnd_imp_model 		= param->getString( "ground_impedence_model" );
	usrattfile 			= param->getString( "use_attn_file" );
	z_min 				= param->getFloat( "zground_km" ) * 1000.0;    // meters
	//freq 				= param->getFloat( "freq" );
	maxrange 			= param->getFloat( "maxrange_km" ) * 1000.0;
	maxheight 			= param->getFloat( "maxheight_km" ) * 1000.0;      // @todo fix elsewhere that m is required
	sourceheight 		= param->getFloat( "sourceheight_km" ) * 1000.0;
	receiverheight 		= param->getFloat( "receiverheight_km" ) * 1000.0;
	tol 				= 1.0e-8;
	Nz_grid 			= param->getInteger( "Nz_grid" );
	Nrng_steps 			= param->getInteger( "Nrng_steps" );
	Lamb_wave_BC 		= param->getBool( "Lamb_wave_BC" );
	write_2D_TLoss  	= param->getBool( "write_2d_tloss" );
	write_phase_speeds 	= param->getBool( "write_phase_speeds" );
  write_atmosphere    = param->getBool( "write_atm_profile" );
	write_modes 		= param->getBool( "write_modes" );
  write_lossless    = param->getBool( "write_lossless" );
	//write_dispersion 	= param->getBool( "write_dispersion" );
  dispersion_file = param->getString( "dispersion_file" );
  append_dispersion_file = param->getBool( "append_dispersion_file" );
	Nby2Dprop 			= param->getBool( "multiprop" );
	turnoff_WKB 		= param->getBool( "turnoff_WKB" );
	z_min_specified     = param->wasFound( "zground_km" );
  user_tag        = param->getString("filetag");
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


	if (write_phase_speeds || write_2D_TLoss || write_modes || (!dispersion_file.empty())) {
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
//    freq = param->getFloat( "freq" );
    f_vec = new double[ 1 ];
    f_vec[ 0 ] = param->getFloat( "freq" );
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
    atm_profile->copy_vector_property( "C0", "_C0_" );
    atm_profile->convert_property_units( "_C0_", Units::fromString( "m/s" ) );

    // do we need to calculate temperature?
    if (!atm_profile->contains_vector("T" )) {
      std::cout << "No temperature provided, calculating from sound speed"
            << std::endl;
      atm_profile->calculate_temperature_from_sound_speed( "T", "_C0_",
        Units::fromString( "K" ) );
    }
  } else {
    if (atm_profile->contains_vector( "P" ) && atm_profile->contains_vector("RHO")) {
      atm_profile->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO", 
        Units::fromString( "m/s" ) );
    } else {
      atm_profile->calculate_sound_speed_from_temperature( "_C0_", "T", 
        Units::fromString( "m/s" ) );
    }
  }
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
  
  atm_profile->calculate_wind_speed( "_WS_", "U", "V" );
	atm_profile->calculate_wind_direction( "_WD_", "U", "V" );
//	if (usrattfile.empty()) {
//    atm_profile->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", freq );
//  } else {
//    atm_profile->read_attenuation_from_file( "_ALPHA_", usrattfile );
//  }

	for (int i=0; i<Nz_grid; i++) {
		Hgt[i] = z_min + i*dz; // Hgt[0] = zground MSL
		rho[i] = atm_profile->get( "RHO", Hgt[i]);
		Pr[i]  = atm_profile->get( "P", Hgt[i] );
		T[i]   = atm_profile->get( "T", Hgt[i] );
		zw[i]  = atm_profile->get( "U", Hgt[i] );
		mw[i]  = atm_profile->get( "V", Hgt[i] );
//		alpha[i]   = atm_profile->get( "_ALPHA_", Hgt[i] );
	}
}

// utility to print the parameters to the screen
void NCPA::WModeSolver::printParams() {
	printf(" Effective Sound Speed Normal Modes Solver Parameters:\n");
	std::ostringstream fstr;
	fstr << "[ ";
	for (int i = 0; i < Nfreq; i++) {
		if (i > 0) {
			fstr << ", ";
		}
		fstr << f_vec[i];
	}
	fstr << " ]";
	printf("                   freq : %s\n", fstr.str().c_str() );
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
  printf("    write_lossless flag : %d\n", write_lossless);
	printf("    write_2d_tloss flag : %d\n", write_2D_TLoss);
	printf("write_phase_speeds flag : %d\n", write_phase_speeds);
	//printf("  write_dispersion flag : %d\n", write_dispersion);
	printf("       write_modes flag : %d\n", write_modes);
	printf("         multiprop flag : %d\n", Nby2Dprop);
	printf("       turnoff_WKB flag : %d\n", turnoff_WKB);
	printf("    atmospheric profile : %s\n", atmosfile.c_str());
  if (!dispersion_file.empty()) {
    printf("        Dispersion file : %s\n", dispersion_file.c_str());
    printf(" Dispersion file append : %d\n", append_dispersion_file );
  }
	if (!usrattfile.empty()) {
		printf("  User attenuation file : %s\n", usrattfile.c_str());
	}


	printf("       wvnum_filter_flg : %d\n", wvnum_filter_flg);
	if (wvnum_filter_flg==1) {
		printf("                  c_min : %g m/s\n", c_min);
		printf("                  c_max : %g m/s\n", c_max);
	}
}


int NCPA::WModeSolver::solve() {
  PetscErrorCode ierr;
  PetscInt nconv;

  int    i, select_modes, nev, it, fi;
  double dz, admittance, h2, rng_step, z_min_km;
  double k_min, k_max;			
  double *diag, *kd, *md, *cd, *kH, *k_s, **v, **v_s;	
  complex<double> *k_pert;

  //alpha  = new double [Nz_grid];
  diag   = new double [Nz_grid];
  kd     = new double [Nz_grid];
  md     = new double [Nz_grid];
  cd     = new double [Nz_grid];
  kH     = new double [MAX_MODES];
  k_s    = new double [MAX_MODES];  
  k_pert = new complex<double> [MAX_MODES];
  v      = dmatrix(Nz_grid,MAX_MODES);
  v_s    = dmatrix(Nz_grid,MAX_MODES); 	

  nev   = 0;
  k_min = 0; 
  k_max = 0;

  rng_step = maxrange/Nrng_steps;  		// range step [meters]
  dz       = (maxheight - z_min)/Nz_grid;	// the z-grid spacing
  h2       = dz*dz;
  //dz_km    = dz/1000.0;
  z_min_km = z_min/1000.0;

  // Initialize Slepc
  SlepcInitialize(PETSC_NULL,PETSC_NULL,(char*)0,PETSC_NULL);
  
  for (fi = 0; fi < Nfreq; fi++) {
    double freq = f_vec[ fi ];
	std::cout << "Calculating at " << freq << " Hz:" << std::endl;
	this->calculateAttenuation( freq );

    //
    // loop over azimuths (if not (N by 2D) it's only one azimuth
    //
    for (it=0; it<Naz; it++) {
      
      azi = azi_min + it*azi_step; // degrees (not radians)
      cout << endl << "Now processing azimuth = " << azi << " (" << it+1 << " of " 
           << Naz << ")" << endl;

      //atm_profile->setPropagationAzimuth(azi);
      atm_profile->calculate_wind_component("_WC_", "_WS_", "_WD_", azi );
      atm_profile->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );
      atm_profile->get_property_vector( "_CEFF_", c_eff );
      atm_profile->add_property( "_AZ_", azi, NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );

      // compute absorption
      //getAbsorption(Nz_grid, dz, atm_profile, freq, usrattfile, alpha);

      //
      // ground impedance model
      //
      // at the ground the BC is: Phi' = (a - 1/2*dln(rho)/dz)*Phi
      // for a rigid ground a=0; and the BC is the Lamb Wave BC:
      // admittance = -1/2*dln(rho)/dz
      //  
      /*
      admittance = 0.0; // default
      if ((gnd_imp_model.compare("rigid")==0) && Lamb_wave_BC) { //Lamb_wave_BC
          admittance = -atm_profile->drhodz(z_min/1000.0)/1000.0/atm_profile->rho(z_min/1000.0)/2.0; // SI units
      }
      else if (gnd_imp_model.compare("rigid")==0) {
          admittance = 0.0; // no Lamb_wave_BC
      }
      else {
          std::ostringstream es;
          es << "This ground impedance model is not implemented yet: " << gnd_imp_model;
          throw invalid_argument(es.str());
      }
      */
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
        admittance, freq, diag, kd, md, cd, &k_min, &k_max, turnoff_WKB);

      // if wavenumber filtering is on, redefine k_min, k_max
      if (wvnum_filter_flg) {
          k_min = 2*PI*freq/c_max;
          k_max = 2*PI*freq/c_min;
      }

      i = getNumberOfModes(Nz_grid,dz,diag,k_min,k_max,&nev);
      
      //
      // set up parameters for eigenvalue estimation; double dimension of problem for linearization
      //
      std::ostringstream oss;
      oss << "______________________________________________________________________" 
          << std::endl << std::endl
          << " -> Solving wide-angle problem at " << freq << " Hz and " << azi << "deg ("
          << nev * 2 << " modes)..." << std::endl
          << " -> Discrete spectrum: " << 2*PI*freq/k_max << " to " << 2*PI*freq/k_min 
          << " m/s" << std::endl
          << " -> Quadratic eigenvalue problem  - double dimensionality.";

      i = NCPA::EigenEngine::doWideAngleCalculation( Nz_grid, dz, k_min, k_max,
        tol, nev, kd, md, cd, &nconv, kH, v, oss.str() );

      // select modes and do perturbation
      doSelect(Nz_grid, nconv, k_min, k_max, kH, v, k_s, v_s, &select_modes);  
      doPerturb(Nz_grid, z_min, dz, select_modes, freq, k_s, v_s, alpha, k_pert);
      //cout << "Found " << select_modes << " relevant modes" << endl;

      //
      // Output data  
      //
      std::string fname_lossy = "", fname_lossless = "";
      if (Nby2Dprop) { // if (N by 2D is requested)
          fname_lossy = tag_filename(NCPAPROP_WMOD_FILENAME_1D_MULTIPROP);
          fname_lossless = "";
          if (write_lossless) {
            fname_lossless = tag_filename(NCPAPROP_WMOD_FILENAME_1D_LOSSLESS_MULTIPROP);
          }
          getTLoss1DNx2(azi, select_modes, dz, Nrng_steps, rng_step, sourceheight, 
          	receiverheight, rho, k_pert, v_s, Nby2Dprop, it,
          	fname_lossy, fname_lossless );
      } else {
          cout << "Writing to file: 1D transmission loss" << endl;
          fname_lossy = tag_filename(NCPAPROP_WMOD_FILENAME_1D);
          std::cout << "Lossy: " << fname_lossy << std::endl;
          fname_lossless = "";
          if (write_lossless) {
            fname_lossless = tag_filename(NCPAPROP_WMOD_FILENAME_1D_LOSSLESS);
            std::cout << "Lossless:" << fname_lossless << std::endl;
          }
          getTLoss1D(select_modes, dz, Nrng_steps, rng_step, sourceheight, receiverheight, 
          	rho, k_pert, v_s, fname_lossy, fname_lossless );

          if (write_2D_TLoss) {
              cout << "Writing to file: 2D transmission loss" << std::endl;
              fname_lossy = tag_filename(NCPAPROP_WMOD_FILENAME_2D);
              std::cout << "Lossy: " << fname_lossy << std::endl;
              fname_lossless = "";
              if (write_lossless) {
                fname_lossless = tag_filename(NCPAPROP_WMOD_FILENAME_2D_LOSSLESS);
                std::cout << "Lossless:" << fname_lossless << std::endl;
              }
              getTLoss2D(Nz_grid,select_modes,dz,Nrng_steps,rng_step,sourceheight,rho,
              	k_pert,v_s, fname_lossy, fname_lossless );
          }

          if (write_phase_speeds) {
              cout << "Writing to file: phase speeds" << std::endl;
              writePhaseSpeeds(select_modes,freq,k_pert,
                tag_filename(NCPAPROP_WMOD_FILENAME_PHASE_SPEEDS));
          }

          if (write_modes) {
              cout << "Writing to file: the modes\n";
              writeEigenFunctions(Nz_grid,select_modes,dz,v_s,"wnm");
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

          if (write_atmosphere && (fi == 0)) {
            std::string atm_output_filename =
              tag_filename(NCPAPROP_WMOD_FILENAME_ATMOSPHERE);
            std::cout << "Writing atmosphere to "
              << atm_output_filename << std::endl;
            std::vector<std::string> keylist;
            keylist.push_back( "U" );
            keylist.push_back( "V" );
            keylist.push_back( "T" );
            keylist.push_back( "RHO" );
            keylist.push_back( "P" );
            keylist.push_back( "_C0_" );
            keylist.push_back( "_CEFF_" );
            std::ofstream atmout( atm_output_filename );
            atm_profile->print_atmosphere( keylist, "Z", atmout );
            atmout.close();
          }
      }
      
      atm_profile->remove_property( "_WC_" );
      atm_profile->remove_property( "_CEFF_" );
      atm_profile->remove_property( "_AZ_" );
    
    } // end loop by azimuths  

  } // end loop by frequencies
  
  // Finalize Slepc
  ierr = SlepcFinalize();CHKERRQ(ierr);
    
  // free the rest of locally dynamically allocated space
  delete[] diag;
  delete[] kd;
  delete[] md;
  delete[] cd;
  delete[] kH;
  delete[] k_s;
  delete[] k_pert;
  free_dmatrix(v, Nz_grid, MAX_MODES);
  free_dmatrix(v_s, Nz_grid, MAX_MODES);

  return 0;
}


int NCPA::WModeSolver::doSelect(int nz, int n_modes, double k_min, double k_max, 
	double *kH, double **v, double *k_s, double **v_s, int *select_modes)
{
  int cnt = -1;
  int i, j;

  for (j=0; j<n_modes; j++) {
      if ((kH[j] >= k_min) && (kH[j] <= k_max)) {
          cnt++;
          for (i=0; i<nz; i++) {
              v_s[i][cnt] = v[i][j];
          }
          k_s[cnt] = kH[j];
      }
  }
  *select_modes = cnt + 1; //+1 because cnt started from zero
  return 0;
}



int NCPA::WModeSolver::getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, 
  double dz, Atmosphere1D *p, double admittance, double freq, double *diag, double *kd, double *md, 
  double *cd, double *k_min, double *k_max, bool turnoff_WKB)
{
  //     use atmospherics input for the trace of the matrix.
  //     the vector diag can be used to solve the modal problem
  //     also returns the bounds on the wavenumber spectrum, [k_min,k_max]
  int i, top;
  double azi_rad, gamma, omega, bnd_cnd, windz, ceffmin, ceffmax, ceff_grnd, cefftop, cz;
  double z_min_km, z_km, dz_km;
  double kk, dkk, k_eff, k_gnd, k_max_full, wkbIntegral, wkbTerm;
  //double rho_factor; 
  double *ceffz;
  ceffz = new double [nz];

  //FILE *profile= fopen("profile.int", "w");

  gamma    = 1.4;
  z_min_km = z_min/1000.0;
  dz_km    = dz/1000.0;
  omega    = 2*PI*freq;

  azi_rad  = NCPA::Units::convert( p->get( "_AZ_" ), NCPA::UNITS_ANGLE_DEGREES, UNITS_ANGLE_RADIANS );

  z_km      = z_min_km; 
  cz        = p->get( "_C0_", z_min );
  windz     = p->get( "_WC_", z_min );
  ceff_grnd = p->get( "_CEFF_", z_min );
  ceffmin   = ceff_grnd;  // in m/s; initialize ceffmin
  ceffmax   = ceffmin;    // initialize ceffmax   

  for (i=0; i<nz; i++) {
     cz         = p->get( "_C0_", Hgt[ i ] );
	    ceffz[ i ] = p->get( "_CEFF_", Hgt[ i ] );
      windz      = p->get( "_WC_", Hgt[ i ] );

	    // we neglect rho_factor - it does not make sense to have in this approximation
	    //rho_factor is 1/2*rho_0"/rho_0 - 3/4*(rho_0')^2/rho_0^2
	    //rho_factor = (-3/4*pow(p->drhodz(z_km)/p->rho(z_km),2) + 1/2*p->ddrhodzdz(z_km)/p->rho(z_km))*1e-6; // in meters^(-2)!
	    kd[i]      = pow(omega/cz,2); // + rho_factor;
	    md[i]      = pow(windz/cz,2) - 1;
	    cd[i]      = -2*omega*(windz/pow(cz,2));
	    diag[i]    = pow(omega/ceffz[ i ],2); // + rho_factor;
	    if (ceffz[ i ] < ceffmin)
	     		ceffmin = ceffz[ i ];     // approximation to find minimum sound speed in problem
	    if (ceffz[ i ] > ceffmax)
	     		ceffmax = ceffz[ i ];
	     
	    //z_km += dz_km; 
  }
  //fclose(profile);

  bnd_cnd = (1./(dz*admittance+1))/(pow(dz,2)); // bnd cnd assuming centered fd
  diag[0] = bnd_cnd + diag[0];
  kd[0]   = bnd_cnd + kd[0];


  // use WKB trick for ground to ground propagation.
  // @todo change this to compare with Z0, not 0.0
  if ((fabs(sourceheight)<1.0e-3) && (fabs(receiverheight)<1.0e-3) && (!turnoff_WKB)) {
      //
      // WKB trick for ground to ground propagation. 
      // Cut off lower phasespeed (highest wavenumber) where tunneling 
      // is insignificant (freq. dependent)
      //
      k_max_full = omega/ceffmin;
      k_gnd      = omega/ceff_grnd;
      dkk        = (pow(k_max_full,2)-pow(k_gnd,2))/100.0;

      //
      // dv: when ceffmin is at the ground dkk can be very small but non-zero (rounding errors?)
      // in that case we need to skip the next for loop; otherwise it will be executed
      // with the small increment dkk
      //
      kk = pow(k_gnd,2); //initialization of kk in case the loop is skipped

      if (dkk >1e-10) { // do this loop only if dkk is nonzero (significantly)
		      for (kk = pow(k_gnd,2); kk < pow(k_max_full,2); kk=kk+dkk) {
              i           = 0;
              wkbIntegral = 0.0;
              wkbTerm     = 1.0; 
              z_km = z_min_km; 
              while (wkbTerm > dkk) {
                  k_eff = omega/ceffz[i];
                  wkbTerm = abs(kk - pow(k_eff,2));
                  wkbIntegral = wkbIntegral + dz*sqrt(wkbTerm); // dz should be in meters
                  i++;
                  z_km += dz_km;
              } 
              if (wkbIntegral >= 10.0) {
                  printf("\n WKB fix: new phasevelocity minimum: %6.2f m/s (was %6.2f m/s); \n WKBIntegral= %12.7f at z = %6.2f km\n", omega/sqrt(kk), omega/k_max_full, wkbIntegral, z_km);
                  break;
              }
          }
      }
      
      *k_max = sqrt(kk);   // use this for ground-to-ground 1D Tloss 
											     //(uses WKB trick to include only non-vanishing modes at the ground)
	    // *k_max = k_max_full;
  }
  else { // not ground-to-ground propagation
      *k_max = omega/ceffmin; // same as k_max_full
  } 

  top     = nz - ((int) nz/10);
  z_km    = z_min_km + (top+1)*dz_km;  
  cz      = p->get( "_C0_", Hgt[ top+1 ] );
  windz   = p->get( "_WC_", Hgt[ top+1 ] );
  cefftop = p->get( "_CEFF_", Hgt[ top+1 ] );
  //cz      = sqrt(gamma*Pr[top+1]/rho[top+1]); // in m/s
  //windz   = zw[top+1]*sin(azi_rad) + mw[top+1]*cos(azi_rad); 
  //cefftop = cz + windz;
  *k_min  = omega/cefftop;
  
  delete [] ceffz;

  return 0;
}
