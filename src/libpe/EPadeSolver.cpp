#include "EPadeSolver.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <string>
#include <vector>
#include <deque>
#include <cfloat>
#include <fstream>
#include <stdexcept>
#include <cstdint>
#include <algorithm>

#include "petscksp.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"
#include "gsl/gsl_blas.h"

#include "ncpaprop_common.h"
#include "ncpaprop_atmosphere.h"
#include "ncpaprop_petsc.h"


#ifndef PI
#define PI 3.14159
#endif

#define RHO_B 5000.0

//////////////////////////////////////////////////////////////////////////////
// Setter methods for calculation parameters
//////////////////////////////////////////////////////////////////////////////

void NCPA::EPadeSolver::set_max_range( double r, const std::string &u ) {
	r_max.set( r, u );
	r_max.convert_units( NCPAPROP_EPADE_PE_UNITS_R );
}

void NCPA::EPadeSolver::set_max_height( double z, const std::string &u ) {
	z_max.set( z, u );
	z_max.convert_units( NCPAPROP_EPADE_PE_UNITS_Z );
}

void NCPA::EPadeSolver::set_source_height( double z, const std::string &u ) {
	z_source.set( z, u );
	z_source.convert_units( NCPAPROP_EPADE_PE_UNITS_Z );
}

void NCPA::EPadeSolver::set_receiver_height( double z, const std::string &u ) {
	z_receiver.set( z, u );
	z_receiver.convert_units( NCPAPROP_EPADE_PE_UNITS_Z );
}

void NCPA::EPadeSolver::set_requested_range_steps( size_t n ) {
	nr_requested = n;
}

void NCPA::EPadeSolver::set_frequency( double f ) {
	f_vector.clear();
	f_vector.push_back( f );
}

void NCPA::EPadeSolver::set_pade_order( size_t n ) {
	pade_order = n;
}

void NCPA::EPadeSolver::set_requested_height_step( double dz, const std::string &u ) {
	dz_requested.set( dz, u );
	dz_requested.convert_units( NCPAPROP_EPADE_PE_UNITS_Z );
}

void NCPA::EPadeSolver::set_starter( NCPA::StarterType st, const std::string &fname ) {
	starter_type = st;
	starter_filename_in = fname;
}

void NCPA::EPadeSolver::set_attenuation( AttenuationType att, const std::string &fname ) {
	attenuation_type = att;
	attenuation_filename_in = fname;
}

void NCPA::EPadeSolver::set_azimuth( double a ) {
	az_vector.clear();
	az_vector.push_back(NCPA::normalizeAzimuth(a));
}

void NCPA::EPadeSolver::set_azimuths( double start_az, double end_az, double daz ) {
	while (start_az > end_az) {
		start_az -= 360.0;
	}
	az_vector.clear();
	for (double a = start_az; a <= end_az; a += daz) {
		az_vector.push_back(NCPA::normalizeAzimuth(a));
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Below here is pre-refactor code /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


NCPA::EPadeSolver::EPadeSolver() {
	set_default_values();
}

void NCPA::EPadeSolver::set_default_values() {
	// actual default values
	absorption_layer_mu = 0.01;
	c_underground = 50000000000.0;
	top_layer = true;
	write1d = true;

	// null values otherwise.  Pointers:
	z = NULL; z_abs = NULL; r = NULL; tl = NULL;
	zgi_r = NULL; atm_profile_2d = NULL;

	// doubles
	z_min = 0.0; z_ground = 0.0;
	z_bottom = 0.0;
	top_layer_thickness_m = -1.0;

	// complex
	user_ground_impedence = 0.0;

	// int
	NZ = 0; NR = 0; nzplot = 0;

	// bool
	use_atm_1d = false; use_atm_2d = false; use_atm_toy = false; use_topo = false;
	z_ground_specified = false;
	write2d = false;
	write_starter = false; write_topo = false;
	user_ground_impedence_found = false; write_atmosphere = false;
	pointsource = true; _write_source_function = false;

	// string
	topofile = "";

	// turbulence
	use_turbulence = false;
	random_turbulence = true;
	turbulence_k1 = 0.1;
	turbulence_k2 = 20.0;
	turbulence_size = 20;
	Lt = 100.0;
	temperature_factor = 1.0e-10;
	velocity_factor    = 1.0e-8;
	// turbulence_vec1 = PETSC_NULL;
}

std::string NCPA::EPadeSolver::tag_filename( std::string basename ) {
	return user_tag + basename;
}


NCPA::EPadeSolver::EPadeSolver( NCPA::ParameterSet *param ) {

	set_default_values();
	set_max_range( param->getFloat( "maxrange_km" ), "km" );

	// obtain the parameter values from the user's options
	// @todo add units to input scalar quantities
//	r_max	 			= param->getFloat( "maxrange_km" ) * 1000.0;
//	z_max.set( param->getFloat( "maxheight_km" ), UNITS_DISTANCE_KILOMETERS );
//  	z_max	 			= param->getFloat( "maxheight_km" ) * 1000.0;      // @todo fix elsewhere that m is required
//  	zs			 		= param->getFloat( "sourceheight_km" ) * 1000.0;
//  	zr 					= param->getFloat( "receiverheight_km" ) * 1000.0;
//  	NR_requested 		= param->getInteger( "Nrng_steps" );
//  	freq 				= param->getFloat( "freq" );
//	npade 				= param->getInteger( "npade" );
//	starter 			= param->getString( "starter" );
//	attnfile 			= param->getString( "attnfile" );
//	user_starter_file   = param->getString( "starterfile" );
	topofile			= param->getString( "topofile" );
	user_tag			= param->getString( "filetag" );
	if (user_tag.size() > 0) {
		user_tag += ".";
	}
	top_layer_thickness_m = param->getFloat( "top_layer_thickness_m" );
	linesourcefile 		= param->getString( "linesourcefile" );
	if (linesourcefile.size() > 0) {
		pointsource = false;
	}

	// flags
//	lossless 				= param->wasFound( "lossless" );
	use_atm_1d				= param->wasFound( "atmosfile" );
	use_atm_2d				= param->wasFound( "atmosfile2d" );
	//use_atm_toy			= param->wasFound( "toy" );
	top_layer				= !(param->wasFound( "disable_top_layer" ));
	use_topo				= param->wasFound( "topo" );
	write2d 				= param->wasFound( "write_2d_tloss" );
	write_starter       	= param->wasFound( "write_starter" );
	_write_source_function 	= param->wasFound( "write_source" );
	write_topo 		 		= param->wasFound( "write_topography" );
	write_atmosphere 		= param->wasFound( "write_atm_profile" );

//	// Handle differences based on single vs multiprop
//	double min_az, max_az, step_az;
//	if (multiprop) {
//		if (use_atm_2d) {
//			std::cerr << "Range-dependent 2-D atmosphere incompatible with multiple azimuth propagation"
//					  << std::endl;
//			exit(0);
//		}
//		if (write2d) {
//			std::cout << "Multi-azimuth propagation requested, disabling 2-D output" << std::endl;
//			write2d = false;
//		}
//		if (use_topo) {
//			std::cout << "Multi-azimuth propagation requested, disabling topography flag" << std::endl;
//			use_topo = false;
//		}
//		min_az 			= param->getFloat( "azimuth_start" );
//		max_az 			= param->getFloat( "azimuth_end" );
//		step_az 		= param->getFloat( "azimuth_step" );
//
//		// set up azimuth vector
//		while (min_az > max_az) {
//			// crossover the zero azimuth point
//			min_az -= 360.0;
//		}
//		NAz 			= (int) ((max_az - min_az)/step_az) + 1;
//
//		// initialize output file
//		std::cout << "Initializing file "
//			<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP) << std::endl;
//		std::ofstream truncator( tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP),
//			std::ofstream::out | std::ofstream::trunc );
//		truncator.close();
//
//	} else {
//		NAz 			= 1;
//		min_az 			= param->getFloat( "azimuth" );
//		max_az 			= min_az;
//		step_az 		= 0;
//	}
//	azi = NCPA::zeros<double>( NAz );
//	for (size_t i = 0; i < NAz; i++) {
//		azi[ i ] = min_az + i * step_az;
//	}

//	if (starter == "user") {
//		if (user_starter_file.size() == 0) {
//			throw std::runtime_error( "No user starter file specified!" );
//		}
//	}

	//NCPA::Atmosphere1D *atm_profile_1d;
	if (use_atm_1d) {
		atm_profile_2d = new NCPA::StratifiedAtmosphere2D(
			param->getString( "atmosfile" ),
			param->getString("atmosheaderfile") );
	} else if (use_atm_2d) {
		atm_profile_2d = new NCPA::ProfileSeriesAtmosphere2D(
			param->getString( "atmosfile2d" ),
			param->getString( "atmosheaderfile" ) );
	} else {
		std::cerr << "Unknown atmosphere option selected" << std::endl;
		exit(0);
	}
	atm_profile_2d->convert_range_units( NCPAPROP_EPADE_PE_UNITS_R );
	if (r_max.get() > atm_profile_2d->get_maximum_valid_range() ) {
		atm_profile_2d->set_maximum_valid_range( r_max.get() );
	}

	// altitude units
	atm_profile_2d->convert_altitude_units( NCPAPROP_EPADE_PE_UNITS_Z );

	// Ground height is treated differently depending on whether we're
	// using topography or not
	if (use_topo) {
		// first, do we get topography from a file?
		if ( topofile.size() > 0 ) {
			atm_profile_2d->remove_property("Z0");
			atm_profile_2d->read_elevation_from_file( topofile );
			z_ground_specified = true;
		} else if (param->wasFound("groundheight_km")) {
			z_ground = param->getFloat( "groundheight_km" ) * 1000.0;
			z_ground_specified = true;
			std::cout << "Overriding profile Z0 value with command-line value " << z_ground
			    << " m" << std::endl;
			atm_profile_2d->remove_property("Z0");
			atm_profile_2d->add_property( "Z0", z_ground,
					NCPAPROP_EPADE_PE_UNITS_Z );
			atm_profile_2d->finalize_elevation_from_profiles();
		} else {
			atm_profile_2d->finalize_elevation_from_profiles();
		}
		z_ground = atm_profile_2d->get_interpolated_ground_elevation(0.0);
	} else {
		// constant elevation
		if (param->wasFound("groundheight_km")) {
			z_ground = param->getFloat( "groundheight_km" ) * 1000.0;
			z_ground_specified = true;

			std::cout << "Overriding profile Z0 value with command-line value " << z_ground
			     << " m" << std::endl;
			atm_profile_2d->remove_property("Z0");
			atm_profile_2d->add_property( "Z0", z_ground, NCPAPROP_EPADE_PE_UNITS_Z );
		} else {
			if (!(atm_profile_2d->contains_scalar(0,"Z0"))) {
				z_ground = atm_profile_2d->get_minimum_altitude( 0.0 );
				atm_profile_2d->add_property("Z0",z_ground,
					atm_profile_2d->get_altitude_units(0.0));
			}
		}
		atm_profile_2d->convert_property_units( "Z0", NCPAPROP_EPADE_PE_UNITS_Z );
		z_ground = atm_profile_2d->get( 0.0, "Z0" );
	}

	// z_min = atm_profile_2d->get_minimum_altitude( 0.0 );
	// z_ground = z_min;




//	// set units
//	if (atm_profile_2d->contains_vector(0,"U")) {
//		atm_profile_2d->convert_property_units( "U", NCPAPROP_EPADE_PE_UNITS_U );
//	}
//	if (atm_profile_2d->contains_vector(0,"V")) {
//		atm_profile_2d->convert_property_units( "V", NCPAPROP_EPADE_PE_UNITS_V );
//	}
//	if (atm_profile_2d->contains_vector(0,"T")) {
//		atm_profile_2d->convert_property_units( "T", NCPAPROP_EPADE_PE_UNITS_T );
//	}
//	if (atm_profile_2d->contains_vector(0,"P")) {
//		atm_profile_2d->convert_property_units( "P", NCPAPROP_EPADE_PE_UNITS_P );
//	}
//
//	// need density
//	if (atm_profile_2d->contains_vector(0,"RHO")) {
//		atm_profile_2d->convert_property_units( "RHO", NCPAPROP_EPADE_PE_UNITS_RHO );
//	} else {
//		std::cout << "No density provided, calculating from temperature and pressure" << std::endl;
//		for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm_profile_2d->first_profile();
//			it != atm_profile_2d->last_profile(); ++it) {
//			if ( (*it)->contains_vector("T") && (*it)->contains_vector("P") ) {
//				(*it)->calculate_density_from_temperature_and_pressure(
//					"RHO", "T", "P", NCPAPROP_EPADE_PE_UNITS_RHO );
//			} else {
//				throw std::runtime_error( "No RHO provided, and at least one of T and P is missing." );
//			}
//		}
//	}
	// z_ground = atm_profile_2d->get( 0.0, "Z0" );

	// calculate derived quantities
//	double c0;
//	for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm_profile_2d->first_profile();
//		 it != atm_profile_2d->last_profile(); ++it) {
//		if ( (*it)->contains_vector( "C0" ) ) {
//			(*it)->convert_property_units( "C0", NCPAPROP_EPADE_PE_UNITS_C );
//			(*it)->copy_vector_property( "C0", "_C0_" );
//			c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
//		} else {
//			if ( (*it)->contains_vector("P") && (*it)->contains_vector("RHO") ) {
//				(*it)->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO",
//						NCPAPROP_EPADE_PE_UNITS_C );
//				c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
//			} else if ( (*it)->contains_vector("T") ) {
//				(*it)->calculate_sound_speed_from_temperature( "_C0_", "T",
//						NCPAPROP_EPADE_PE_UNITS_C );
//				c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
//			} else if ( (*it)->contains_vector( "CEFF" ) ) {
//				c0 = atm_profile_2d->get( 0.0, "CEFF", z_ground );
//			} else {
//				throw std::runtime_error( "Cannot calculate static sound speed: None of CEFF, C0, T, or (P and RHO) are specified in atmosphere.");
//			}
//		}
//	}

//	// wind speed
//	if (atm_profile_2d->contains_vector(0,"WS")) {
//		atm_profile_2d->convert_property_units("WS", NCPAPROP_EPADE_PE_UNITS_U );
//		atm_profile_2d->copy_vector_property( "WS", "_WS_" );
//	} else if (atm_profile_2d->contains_vector(0,"U")
//		&& atm_profile_2d->contains_vector(0,"V")) {
//		atm_profile_2d->calculate_wind_speed( "_WS_", "U", "V" );
//	}
//
//	// wind direction
//	if (atm_profile_2d->contains_vector(0,"WD")) {
//		atm_profile_2d->convert_property_units("WD",
//			NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
//		atm_profile_2d->copy_vector_property( "WD", "_WD_" );
//	} else if (atm_profile_2d->contains_vector(0,"U")
//		&& atm_profile_2d->contains_vector(0,"V")) {
//		atm_profile_2d->calculate_wind_direction( "_WD_", "U", "V" );
//	}

	// calculate/check z resolution
//	dz = param->getFloat( "dz_m" );
//	double f0 = f_vector[0];
//	double lambda0 = c0 / f_vector[0];
//  	if (dz <= 0.0) {
//  		dz = lambda0 / 20.0;
//  		double nearestpow10 = std::pow( 10.0, (double)std::floor( (double)std::log10( dz ) ) );
//  		double factor = std::floor( dz / nearestpow10 );
//  		dz = nearestpow10 * factor;
//  		std::cout << "Setting dz to " << dz << " m" << std::endl;
//  	}
//  	if (dz > (c0 / f_vector[0] / 10.0) ) {
//  		std::ostringstream oss;
//		oss << "Altitude resolution is too coarse.  Must be <= " << lambda0 / 10.0 << " meters.";
//  		throw std::runtime_error( oss.str() );
//  	}

  	// calculate ground impedence
  	if (param->wasFound( "ground_impedence_real" ) || param->wasFound( "ground_impedence_imag" ) ) {
  		user_ground_impedence.real( param->getFloat( "ground_impedence_real" ) );
  		user_ground_impedence.imag( param->getFloat( "ground_impedence_imag" ) );
  		user_ground_impedence_found = true;
  	}

  	// create turbulence
	use_turbulence = param->wasFound( "turbulence" );
	turbulence_size = (size_t)(param->getInteger( "n_turbulence" ));
	random_turbulence = !(param->wasFound("turbulence_file"));
	if (!random_turbulence)
		turbulence_file = param->getString("turbulence_file");
	// T0 = param->getFloat( "turbulence_ref_temp" );
	Lt = param->getFloat( "turbulence_scale_m" );
	temperature_factor = param->getFloat( "turbulence_t_factor" );
	velocity_factor    = param->getFloat( "turbulence_v_factor" );
}

NCPA::EPadeSolver::~EPadeSolver() {
	if (atm_profile_2d != NULL)
		delete atm_profile_2d;
}

bool NCPA::EPadeSolver::finalize() {

	// Handle differences based on single vs multiprop
	if (az_vector.size() > 1) {
		if (use_atm_2d) {
			std::cerr << "Range-dependent 2-D atmosphere incompatible with multiple azimuth propagation"
					  << std::endl;
			exit(0);
		}
		if (write2d) {
			std::cout << "Multi-azimuth propagation requested, disabling 2-D output" << std::endl;
			write2d = false;
		}
		if (use_topo) {
			std::cout << "Multi-azimuth propagation requested, disabling topography flag" << std::endl;
			use_topo = false;
		}

		// initialize output file
		std::cout << "Initializing file "
			<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP) << std::endl;
		std::ofstream truncator( tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP),
			std::ofstream::out | std::ofstream::trunc );
		truncator.close();
	}

	// set units
	if (atm_profile_2d->contains_vector(0,"U")) {
		atm_profile_2d->convert_property_units( "U", NCPAPROP_EPADE_PE_UNITS_U );
	}
	if (atm_profile_2d->contains_vector(0,"V")) {
		atm_profile_2d->convert_property_units( "V", NCPAPROP_EPADE_PE_UNITS_V );
	}
	if (atm_profile_2d->contains_vector(0,"T")) {
		atm_profile_2d->convert_property_units( "T", NCPAPROP_EPADE_PE_UNITS_T );
	}
	if (atm_profile_2d->contains_vector(0,"P")) {
		atm_profile_2d->convert_property_units( "P", NCPAPROP_EPADE_PE_UNITS_P );
	}

	// need density
	if (atm_profile_2d->contains_vector(0,"RHO")) {
		atm_profile_2d->convert_property_units( "RHO", NCPAPROP_EPADE_PE_UNITS_RHO );
	} else {
		std::cout << "No density provided, calculating from temperature and pressure" << std::endl;
		for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm_profile_2d->first_profile();
			it != atm_profile_2d->last_profile(); ++it) {
			if ( (*it)->contains_vector("T") && (*it)->contains_vector("P") ) {
				(*it)->calculate_density_from_temperature_and_pressure(
					"RHO", "T", "P", NCPAPROP_EPADE_PE_UNITS_RHO );
			} else {
				throw std::runtime_error( "No RHO provided, and at least one of T and P is missing." );
			}
		}
	}

	// wind speed
	if (atm_profile_2d->contains_vector(0,"WS")) {
		atm_profile_2d->convert_property_units("WS", NCPAPROP_EPADE_PE_UNITS_U );
		atm_profile_2d->copy_vector_property( "WS", "_WS_" );
	} else if (atm_profile_2d->contains_vector(0,"U")
		&& atm_profile_2d->contains_vector(0,"V")) {
		atm_profile_2d->calculate_wind_speed( "_WS_", "U", "V" );
	}

	// wind direction
	if (atm_profile_2d->contains_vector(0,"WD")) {
		atm_profile_2d->convert_property_units("WD",
			NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
		atm_profile_2d->copy_vector_property( "WD", "_WD_" );
	} else if (atm_profile_2d->contains_vector(0,"U")
		&& atm_profile_2d->contains_vector(0,"V")) {
		atm_profile_2d->calculate_wind_direction( "_WD_", "U", "V" );
	}

	// attenuation
	switch (attenuation_type) {
		case NCPA::AttenuationType::USER:
			atm_profile_2d->read_attenuation_from_file( "_ALPHA_", attenuation_filename_in );
			break;
		case NCPA::AttenuationType::SUTHERLAND_BASS:
			if (  !(atm_profile_2d->contains_vector(0.0, "T")
					&& atm_profile_2d->contains_vector(0.0, "P")
					&& atm_profile_2d->contains_vector(0.0, "RHO") ) ) {
				std::cout << "At least one of T, P, or RHO is absent, switching to lossless propagation"
						  << std::endl;
				attenuation_type = NCPA::AttenuationType::NONE;
			}
			break;
		default:
			break;
	}

	if (z_max.get_units() == NCPA::UNITS_NONE) {
		double minmax, maxmax;
		atm_profile_2d->get_maximum_altitude_limits( minmax, maxmax );
		z_max.set( minmax, atm_profile_2d->get_altitude_units(0.0) );
	}

	_ready = true;
	std::ostringstream errors;

	// do checks
	if (atm_profile_2d == NULL) {
		errors << "No valid atmosphere provided!" << std::endl;
		_ready = false;
	}
	if (starter_type == NCPA::StarterType::NONE) {
		errors << "No valid starter indicated!" << std::endl;
		_ready = false;
	}
	if (r_max.get_units() == NCPA::UNITS_NONE) {
		errors << "No propagation range set!" << std::endl;
		_ready = false;
	}
	if (az_vector.size() == 0) {
		errors << "No azimuths set!" << std::endl;
		_ready = false;
	}
	if (!_ready) {
		throw std::runtime_error( errors.str() );
	}
	return _ready;
}

bool NCPA::EPadeSolver::ready() {
	return _ready;
}

int NCPA::EPadeSolver::solve() {
	if (!_ready) {
		throw std::runtime_error("Solver not ready!");
	}
	if (use_topo) {
		return solve_with_topography();
	} else {
		return solve_without_topography();
	}
}

int NCPA::EPadeSolver::solve_without_topography() {
	size_t i;
	std::complex<double> I( 0.0, 1.0 );
	PetscErrorCode ierr;
	PetscInt *indices;
	PetscScalar hank, *contents;
	Mat B, C;   // , q;
	Mat *qpowers = PETSC_NULL, *qpowers_starter = PETSC_NULL;
	Vec psi_o, Bpsi_o; //, psi_temp;
	KSP ksp;

	// for turbulence, if needed
	double *mu_r, *mu_rpdr;
	std::vector<double> rand1, rand2;

	// PC pc;

	// set up z grid for flat ground.  When we add terrain we will need to move this inside
	// the range loop
	//int profile_index = atm_profile_2d->get_profile_index( 0.0 );
	int profile_index;
	double minlimit, maxlimit;
	atm_profile_2d->get_maximum_altitude_limits( minlimit, maxlimit );
//	z_max = NCPA::min( z_max, minlimit );    // lowest valid top value
	z_max.set_value( NCPA::min( z_max.get(), minlimit ) );
	int ground_index = 0;
	std::complex<double> ground_impedence_factor( 0.0, 0.0 );

	// truncate multiprop file if needed
	if (write2d) {
		std::ofstream ofs( tag_filename(NCPAPROP_EPADE_PE_FILENAME_2D),
			std::ofstream::out | std::ofstream::trunc );
		ofs.close();
	}

	// truncate 1-D if necessary
//	if (broadband) {
//		std::ofstream ofs( tag_filename(NCPAPROP_EPADE_PE_FILENAME_1D),
//			std::ofstream::out | std::ofstream::trunc );
//		ofs.close();
//	}

	if (r_max.get() > atm_profile_2d->get_maximum_valid_range() ) {
		atm_profile_2d->set_maximum_valid_range( r_max.get() );
	}

	double c0;
	for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm_profile_2d->first_profile();
		 it != atm_profile_2d->last_profile(); ++it) {
		if ( (*it)->contains_vector( "C0" ) ) {
			(*it)->convert_property_units( "C0", NCPAPROP_EPADE_PE_UNITS_C );
			(*it)->copy_vector_property( "C0", "_C0_" );
			c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
		} else {
			if ( (*it)->contains_vector("P") && (*it)->contains_vector("RHO") ) {
				(*it)->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO",
						NCPAPROP_EPADE_PE_UNITS_C );
				c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
			} else if ( (*it)->contains_vector("T") ) {
				(*it)->calculate_sound_speed_from_temperature( "_C0_", "T",
						NCPAPROP_EPADE_PE_UNITS_C );
				c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
			} else if ( (*it)->contains_vector( "CEFF" ) ) {
				c0 = atm_profile_2d->get( 0.0, "CEFF", z_ground );
			} else {
				throw std::runtime_error( "Cannot calculate static sound speed: None of CEFF, C0, T, or (P and RHO) are specified in atmosphere.");
			}
		}
	}

	// calculate/check z resolution
//	dz = param->getFloat( "dz_m" );
//	double f0 = f_vector[0];
	// @todo move this into f loop for future broadband use
	double dz;
	double lambda0 = c0 / f_vector[0];
	if (dz_requested.get() <= 0.0) {
		dz = lambda0 / 20.0;
		double nearestpow10 = std::pow( 10.0, (double)std::floor( (double)std::log10( dz ) ) );
		double factor = std::floor( dz / nearestpow10 );
		dz = nearestpow10 * factor;
		std::cout << "Setting dz to " << dz << " m" << std::endl;
	}
	if (dz > (c0 / f_vector[0] / 10.0) ) {
		std::ostringstream oss;
		oss << "Altitude resolution is too coarse.  Must be <= " << lambda0 / 10.0 << " meters.";
		throw std::runtime_error( oss.str() );
	}

	atm_profile_2d->get_minimum_altitude_limits( minlimit, z_min );
	//z_min = atm_profile_2d->get_highest_minimum_altitude();
	// if ( (!z_ground_specified) && atm_profile_2d->contains_scalar( 0.0, "Z0" )) {
	// if ( !z_ground_specified ) {
		// z_ground = atm_profile_2d->get( 0.0, "Z0" );
	// 	z_ground = z_min;
	// }
	// z_ground = atm_profile_2d->get(0.0,"Z0");
	if (z_ground < z_min) {
		std::cerr << "Supplied ground height is outside of atmospheric specification." << std::endl;
		exit(0);
	}
  	z_bottom = z_min;
	// fill and convert to SI units
	//double dz       = (z_max - z_ground)/(NZ - 1);	// the z-grid spacing
	NZ = ((int)std::floor((z_max.get() - z_ground) / dz)) + 1;
	z = NCPA::zeros<double>( NZ );
	z_abs = NCPA::zeros<double>( NZ );
	indices = NCPA::zeros<PetscInt>( NZ );
	for (i = 0; i < NZ; i++) {
		z[ i ]     = ((double)i) * dz;
		z_abs[ i ] = z[ i ] + z_ground;
		indices[ i ] = i;
	}
	size_t zr_i = NCPA::find_closest_index<double>( z, NZ, z_receiver.get() );
//	std::cout << "zr_i for zr=" << zr << " is " << zr_i << std::endl;

	if (use_turbulence) {
		if (!random_turbulence) {
			// if (verbose) {
				std::cout << "Reading " << 2*turbulence_size
						  << " values from " << turbulence_file
						  << std::endl;
			// }
			std::ifstream rand_in( turbulence_file );
			if (!rand_in.good()) {
				throw std::runtime_error(
					"Error opening " + turbulence_file);
			}
			rand1.reserve( turbulence_size );
			for (i = 0; i < turbulence_size; i++) {
				rand_in >> rand1[ i ];
				if (!rand_in.good()) {
					std::ostringstream oss;
					oss << "Error reading turbulence numbers from "
						<< turbulence_file;
					throw std::runtime_error(oss.str());
				}
			}
			rand2.reserve( turbulence_size );
			for (i = 0; i < turbulence_size; i++) {
				rand_in >> rand2[ i ];
			}
			rand_in.close();
		}
	}
	
	// constants for now
	double h = dz;
	double h2 = h * h;
	double dr;

	// set up for source atmosphere
	double k0 = 0.0; c0 = 0.0;
	double *c = NCPA::zeros<double>( NZ );
	double *a_t = NCPA::zeros<double>( NZ );
	std::complex<double> *k = NCPA::zeros<std::complex<double>>( NZ );
	std::complex<double> *n = NCPA::zeros<std::complex<double>>( NZ );

	std::complex<double> *source = NCPA::zeros<std::complex<double>>( NZ );
	if (starter_type == NCPA::StarterType::SELF) {
		if (pointsource) {
			std::cout << "Generating point source at " << z_source.get() << "m" << std::endl;
			make_point_source( NZ, z, z_source.get(), 0.0, source );
		} else {
			std::cout << "Reading line source from " << linesourcefile << std::endl;
			read_line_source_from_file( NZ, z, z_ground,
				linesourcefile, source );
		}

		// output source file for checking
		if (_write_source_function) {
			std::cout << "Writing source function to "
					<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_SOURCE) << std::endl;
			write_source(tag_filename(NCPAPROP_EPADE_PE_FILENAME_SOURCE), source, z, NZ);
		}
	}

	// write broadband header for testing
//	if (broadband) {
//		write_broadband_header( tag_filename(NCPAPROP_EPADE_PE_FILENAME_BROADBAND),
//			azi, NAz, f, NF, 1.0e8 );
//	}

	// freq and calc_az hold the current values of azimuth and frequency, respectively
	// these are used in the output routines, so make sure they get set correctly
	// whenever you change frequencies and azimuths
//	for (size_t azind = 0; azind < NAz; azind++) {
	for (std::vector<double>::const_iterator az = az_vector.cbegin(); az != az_vector.cend(); ++az) {
//		calc_az = azi[ azind ];

		profile_index = -1;
		calculate_effective_sound_speed( atm_profile_2d, *az, "_CEFF_" );
		// atm_profile_2d->calculate_wind_component( "_WC_", "_WS_", "_WD_",
		// 	calc_az );
		// atm_profile_2d->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );

		for (std::vector<double>::const_iterator freq = f_vector.cbegin();
				freq != f_vector.cend(); ++freq) {

//			freq = f[ freqind ];
			std::cout << "Infrasound PE code at f = " << *freq << " Hz, azi = "
						<< *az << " deg" << std::endl;

//			if ( (!lossless) && (attnfile.length() == 0)) {
			if (attenuation_type == NCPA::AttenuationType::SUTHERLAND_BASS) {
				atm_profile_2d->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", *freq );
			}

			if (nr_requested == 0) {
		  		dr = 340.0 / *freq;
				NR = (int)std::ceil( r_max.get() / dr );
		  	} else {
		  		NR = nr_requested;
		  		dr = r_max.get() / NR;
		  	}

		  	std::cout << "Setting dr to " << dr << " meters." << std::endl;

		  	r = NCPA::zeros<double>( NR );
		  	zgi_r = NCPA::zeros<int>( NR );

		  	for (i = 0; i < NR; i++) {
		  		r[ i ] = ((double)(i+1)) * dr;
		  	}		
			tl = NCPA::cmatrix( NZ, NR-1 );

			// calculate ground impedence (script A in notes in eq. 12)
			double rho0 = atm_profile_2d->get( 0.0, "RHO", z_ground );
			double lambBC = atm_profile_2d->get_first_derivative( 0.0, "RHO", z_ground ) / (2.0 * rho0);
			//lambBC = 0.0;
			if (user_ground_impedence_found) {
				ground_impedence_factor = *freq * I * 2.0 * PI * rho0 / user_ground_impedence + lambBC;
				std::cout << "Using user ground impedence of " << user_ground_impedence << std::endl;
				//		<< " results in calculated A factor of " << ground_impedence_factor << std::endl;
			} else {
				ground_impedence_factor.real( lambBC );
				ground_impedence_factor.imag( 0.0 );
				std::cout << "Using default rigid ground with Lamb BC" << std::endl;
				//: A factor = " << ground_impedence_factor << std::endl;
			}

			//std::cout << "Using atmosphere index " << profile_index << std::endl;
			calculate_atmosphere_parameters( atm_profile_2d, NZ, z, 0.0, z_ground, attenuation_type,
				top_layer, *freq, use_topo, k0, c0, c, a_t, k, n );

			// calculate turbulence
			if (use_turbulence) {
				mu_r = NCPA::zeros<double>( NZ );
				mu_rpdr = NCPA::zeros<double>( NZ );
				setup_turbulence(rand1, rand2);
			}

			// calculate q matrices
			Mat q;
			build_operator_matrix_without_topography( NZ, z, k0, h2, ground_impedence_factor,
				n, pade_order+1, 0, &q );
			// outputSparseMat( q, NZ, "Q_orig.dat" );
			create_matrix_polynomial( pade_order+1, &q, &qpowers );
			ierr = MatDestroy( &q );CHKERRQ(ierr);

			if (starter_type == NCPA::StarterType::SELF) {
				Mat q_starter;
				build_operator_matrix_without_topography( NZ, z, k0, h2, ground_impedence_factor,
					n, pade_order+1, 0, &q_starter );
				create_matrix_polynomial( pade_order+1, &q_starter, &qpowers_starter );

				// source function calculated earlier
				get_starter_self( NZ, z, source, k0, qpowers_starter,
					pade_order, &psi_o );
				
				ierr = MatDestroy( &q_starter );CHKERRQ(ierr);
			} else if (starter_type == NCPA::StarterType::GAUSSIAN) {
				qpowers_starter = qpowers;
				get_starter_gaussian( NZ, z, z_source.get(), k0, ground_index, &psi_o );
			} else if (starter_type == NCPA::StarterType::USER) {
				get_starter_user( starter_filename_in, NZ, z, &psi_o );
			} else {
				std::cerr << "Unrecognized starter type" << std::endl;
				exit(0);
			}

			if (write_starter) {
				std::cout << "Outputting starter..." << std::endl;
				NCPA::outputVec( psi_o, z, NZ, tag_filename(NCPAPROP_EPADE_PE_FILENAME_STARTER) );
			}

			std::cout << "Finding ePade coefficients..." << std::endl;
			std::vector< std::complex<double> > P, Q;
			std::vector< PetscScalar > taylor = taylor_exp_id_sqrt_1pQ_m1( 2*pade_order, k0*dr );
			calculate_pade_coefficients( &taylor, pade_order, pade_order+1, &P, &Q );
			generate_polymatrices( qpowers_starter, pade_order, NZ, P, Q, &B, &C );

			std::cout << "Marching out field..." << std::endl;
			ierr = VecDuplicate( psi_o, &Bpsi_o );CHKERRQ(ierr);
			contents = NCPA::zeros<PetscScalar>( NZ );

			ierr = KSPCreate( PETSC_COMM_SELF, &ksp );CHKERRQ(ierr);
			ierr = KSPSetOperators( ksp, C, C );CHKERRQ(ierr);
			ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
			for (size_t ir = 0; ir < (NR-1); ir++) {

				double rr = r[ ir ];
				// check for atmosphere change
				if (((int)(atm_profile_2d->get_profile_index( rr ))) != profile_index) {
				
					profile_index = atm_profile_2d->get_profile_index( rr );
					calculate_atmosphere_parameters( atm_profile_2d, NZ, z, rr, z_ground, attenuation_type, top_layer, *freq,
						use_topo, k0, c0, c, a_t, k, n );
					delete_matrix_polynomial( pade_order+1, &qpowers );
					build_operator_matrix_without_topography( NZ, z, k0, h2, 
						ground_impedence_factor, n, pade_order+1, 0, &q );
					create_matrix_polynomial( pade_order+1, &q, &qpowers );
					ierr = MatDestroy( &q );

					taylor.clear();
					taylor = taylor_exp_id_sqrt_1pQ_m1( 2*pade_order, k0*dr );
					calculate_pade_coefficients( &taylor, pade_order, pade_order+1, &P, &Q );
					ierr = MatZeroEntries( B );CHKERRQ(ierr);
					ierr = MatZeroEntries( C );CHKERRQ(ierr);
					generate_polymatrices( qpowers, pade_order, NZ, P, Q, &B, &C );
					std::cout << "Switching to atmosphere index " << profile_index 
						<< " at range = " << rr/1000.0 << " km" << std::endl;
				}

				// get values for current step
				ierr = VecGetValues( psi_o, NZ, indices, contents );CHKERRQ(ierr);

				// apply turbulence
				if (use_turbulence) {
					if (ir == 0) {
						// calculate first step
						calculate_turbulence( rr, NZ, z, k0, 0, mu_r );
					} else {
						std::memcpy( mu_r, mu_rpdr, NZ*sizeof(double) );
					}
					calculate_turbulence( rr + dr, NZ, z, k0, 0, mu_rpdr );

					// apply the turbulent fluctuations.  Do this inside
					// the if() because we need to keep these modifications
					// to psi_o, as opposed to the scaling by the Hankel
					// function below
					for (i = 0; i < NZ; i++) {
						contents[ i ] *= std::exp( I * k0 * dr * 0.5 *
							(mu_r[ i ] + mu_rpdr[ i ]) );
					}

					// store the modified field
					ierr = VecSetValues( psi_o, NZ, indices, contents,
						INSERT_VALUES );CHKERRQ(ierr);
					ierr = VecAssemblyBegin( psi_o );CHKERRQ(ierr);
					ierr = VecAssemblyEnd( psi_o );CHKERRQ(ierr);
				}

				hank = std::sqrt( 2.0 / ( PI * k0 * rr ) ) * exp( I * ( k0 * rr - PI/4.0 ) );
				for (i = 0; i < NZ; i++) {
					tl[ i ][ ir ] = contents[ i ] * hank;
				}
				zgi_r[ ir ] = zr_i;  // constant receiver height

				if ( std::fmod( rr, 1.0e5 ) < dr) {
					std::cout << " -> Range " << rr/1000.0 << " km" << std::endl;
				}

				ierr = MatMult( B, psi_o, Bpsi_o );CHKERRQ(ierr);
				ierr = KSPSetOperators( ksp, C, C );CHKERRQ(ierr);  // may not be necessary
				ierr = KSPSolve( ksp, Bpsi_o, psi_o );CHKERRQ(ierr);
			}
			std::cout << "Stopped at range " << r[ NR-1 ]/1000.0 << " km" << std::endl;

			

			if (az_vector.size() > 1) {
				if (write1d) {
					std::cout << "Writing 1-D output to "
						<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP) << std::endl;
					output1DTL( tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP), *az, true );
				}
			} else { 
				if (write1d) {
					std::cout << "Writing 1-D output to "
						<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_1D) << std::endl;
					output1DTL( tag_filename(NCPAPROP_EPADE_PE_FILENAME_1D), *az, false );
				}
				if (write2d) {
					std::cout << "Writing 2-D output to "
						<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_2D) << std::endl;
					output2DTL( tag_filename(NCPAPROP_EPADE_PE_FILENAME_2D) );
				}
			}

			// write broadband body for testing
//			if (broadband) {
//				write_broadband_results(
//					tag_filename(NCPAPROP_EPADE_PE_FILENAME_BROADBAND),
//					calc_az, freq, r, NR, z_abs, NZ, tl, 1.0e8 );
//			}

			if (write_atmosphere) {
				std::cout << "Writing source atmosphere to "
						<< tag_filename("atm_profile.pe") << std::endl;
				std::vector<std::string> keylist;
				keylist.push_back( "U" );
				keylist.push_back( "V" );
				keylist.push_back( "T" );
				keylist.push_back( "RHO" );
				keylist.push_back( "P" );
				keylist.push_back( "_C0_" );
				keylist.push_back( "_CEFF_" );
				std::ofstream atmout( tag_filename( "atm_profile.pe" ) );
				atm_profile_2d->print_atmosphere( keylist, 0.0, "Z", atmout );
				atmout.close();
			}
			
			std::cout << std::endl;

			// clean up
			if (use_turbulence) {
				delete [] mu_r;
				delete [] mu_rpdr;
				cleanup_turbulence();
			}

			delete_matrix_polynomial( pade_order+1, &qpowers );

			if (starter_type == NCPA::StarterType::SELF) {
				delete_matrix_polynomial( pade_order+1, &qpowers_starter );
				// delete [] qpowers_starter;
			}
//			if (attnfile.length() == 0) {
			if (attenuation_type != NCPA::AttenuationType::USER) {
				atm_profile_2d->remove_property( "_ALPHA_" );
			}
			delete [] r;
			delete [] zgi_r;
			NCPA::free_cmatrix( tl, NZ, NR-1 );
		}
		
		atm_profile_2d->remove_property( "_CEFF_" );
		// atm_profile_2d->remove_property( "_WC_" );
	}

	ierr = MatDestroy( &B );       CHKERRQ(ierr);
	ierr = MatDestroy( &C );       CHKERRQ(ierr);
	ierr = VecDestroy( &psi_o );   CHKERRQ(ierr);
	ierr = VecDestroy( &Bpsi_o );  CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp );     CHKERRQ(ierr);
	
	delete [] k;
	delete [] n;
	delete [] c;
	delete [] a_t;
	delete [] contents;
	delete [] indices;
	delete [] z;
	delete [] z_abs;

	return 1;
}

int NCPA::EPadeSolver::get_starter_user( std::string filename, int NZ, double *z, 
										 Vec *psi ) {

	std::ifstream in( filename );
	std::string line;
	std::vector< std::string > filelines, fields;
//	std::deque< double > z_file, r_file, i_file;
	size_t nlines, i;
	std::ostringstream oss;
	PetscErrorCode ierr;

	std::getline( in, line );

	while ( in.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = NCPA::deblank( line );
		if (line[ 0 ] != '#') {
			filelines.push_back( line );
		}

		std::getline( in, line );
	}
	in.close();

	nlines = filelines.size();
//	double this_z, this_r, this_i;
	double *this_z = NCPA::zeros<double>( nlines );
	std::complex<double> *this_c = NCPA::zeros<std::complex<double>>( nlines );

	for (i = 0; i < nlines; i++) {
		fields = NCPA::split( NCPA::deblank( filelines[ i ] ), " ," );
		if (fields.size() != 3) {
			oss << "EPadeSolver - Error parsing starter line:" << std::endl << line << std::endl 
				<< "Must be formatted as:" << std::endl
				<< "altitude  realpart imagpart" << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		try {
			this_z[ i ] = std::stod( fields[ 0 ] ) * 1000.0;
			this_c[ i ].real( std::stod( fields[ 1 ] ) );
			this_c[ i ].imag( std::stod( fields[ 2 ] ) );
		} catch ( std::invalid_argument &e ) {
			oss << "EPadeSolver - Error parsing starter line:" << std::endl << line << std::endl 
				<< "All fields must be numeric" << std::endl;
			throw std::invalid_argument( oss.str() );
		}
//		z_file.push_back( this_z * 1000.0 );
//		r_file.push_back( this_r );
//		i_file.push_back( this_i );
	}

	ierr = VecCreate( PETSC_COMM_SELF, psi );CHKERRQ(ierr);
	ierr = VecSetSizes( *psi, PETSC_DECIDE, NZ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( *psi ); CHKERRQ(ierr);
	ierr = VecSet( *psi, 0.0 );
	PetscInt *z_indices = NCPA::index_vector<PetscInt>( NZ );
	PetscScalar *new_c = NCPA::zeros<PetscScalar>( NZ );
	interpolate_complex( nlines, this_z, this_c, NZ, z, new_c );
	ierr = VecSetValues( *psi, NZ, z_indices, new_c, INSERT_VALUES );CHKERRQ(ierr);
	ierr = VecAssemblyBegin( *psi );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( *psi );CHKERRQ(ierr);

	delete [] this_z;
	delete [] this_c;
	return 1;
}

void NCPA::EPadeSolver::interpolate_complex( size_t NZ_orig,
		double *z_orig, std::complex<double> *c_orig,
		size_t NZ_new, double *z_new, std::complex<double> *c_new ) {
	double *r_orig = NCPA::zeros<double>( NZ_orig );
	double *i_orig = NCPA::zeros<double>( NZ_orig );
	for (size_t i = 0; i < NZ_orig; i++) {
		r_orig[ i ] = c_orig[ i ].real();
		i_orig[ i ] = c_orig[ i ].imag();
	}
	interpolate_complex( NZ_orig, z_orig, r_orig, i_orig, NZ_new, z_new, c_new );
	delete [] r_orig;
	delete [] i_orig;
}

void NCPA::EPadeSolver::interpolate_complex( size_t NZ_orig,
		double *z_orig, double *r_orig, double *i_orig,
		size_t NZ_new, double *z_new, std::complex<double> *c_new ) {

	gsl_interp_accel *accel_r_, *accel_i_;
	gsl_spline *spline_r_, *spline_i_;
	std::complex<double> J( 0.0, 1.0 );

	accel_r_ = gsl_interp_accel_alloc();
#if GSL_MAJOR_VERSION > 1
	spline_r_ = gsl_spline_alloc( gsl_interp_steffen, NZ_orig );
#else
	spline_r_ = gsl_spline_alloc( gsl_interp_cspline, NZ_orig );
#endif
	gsl_spline_init( spline_r_, z_orig, r_orig, NZ_orig );
	accel_i_ = gsl_interp_accel_alloc();
#if GSL_MAJOR_VERSION > 1
	spline_i_ = gsl_spline_alloc( gsl_interp_steffen, NZ_orig );
#else
	spline_i_ = gsl_spline_alloc( gsl_interp_cspline, NZ_orig );
#endif
	gsl_spline_init( spline_i_, z_orig, i_orig, NZ_orig );

	for (size_t i = 0; i < NZ_new; i++) {
		if (z_new[ i ] < z_orig[ 0 ] || z_new[ i ] > z_orig[ NZ_orig-1 ]) {
			c_new[ i ] = 0.0;
		} else {
			c_new[ i ] = gsl_spline_eval( spline_r_, z_new[ i ], accel_r_ )
				+ J * gsl_spline_eval( spline_i_, z_new[ i ], accel_i_ );
		}
	}

	gsl_spline_free( spline_r_ );
	gsl_spline_free( spline_i_ );
	gsl_interp_accel_free( accel_r_ );
	gsl_interp_accel_free( accel_i_ );
}

//int NCPA::EPadeSolver::interpolate_starter(
//			std::deque<double> &z_orig, std::deque<double> &r_orig,
//			std::deque<double> &i_orig,
//			size_t NZ_new, double *z_new, Vec *psi ) {
//
//	PetscScalar J( 0.0, 1.0 );
//	PetscErrorCode ierr;
//	gsl_interp_accel *accel_r_, *accel_i_;
//	gsl_spline *spline_r_, *spline_i_;
//	double *z_spline, *r_spline, *i_spline;
//
//	PetscInt ii;
//	PetscScalar tempval;
//
//
//	double dz = z_orig[ 1 ] - z_orig[ 0 ];
//	while (z_orig.front() > z_new[ 0 ]) {
//		z_orig.push_front( z_orig[ 0 ] - dz );
//		r_orig.push_front( 0.0 );
//		i_orig.push_front( 0.0 );
//	}
//	while (z_orig.back() < z_new[ NZ_new-1 ]) {
//		z_orig.push_back( z_orig.back() + dz );
//		r_orig.push_back( 0.0 );
//		i_orig.push_back( 0.0 );
//	}
//
//	z_spline = NCPA::zeros<double>( z_orig.size() );
//	r_spline = NCPA::zeros<double>( z_orig.size() );
//	i_spline = NCPA::zeros<double>( z_orig.size() );
//	for (ii = 0; ii < (PetscInt)(z_orig.size()); ii++) {
//		z_spline[ ii ] = z_orig[ ii ];
//		r_spline[ ii ] = r_orig[ ii ];
//		i_spline[ ii ] = i_orig[ ii ];
//	}
//
//	accel_r_ = gsl_interp_accel_alloc();
//#if GSL_MAJOR_VERSION > 1
//	spline_r_ = gsl_spline_alloc( gsl_interp_steffen, z_orig.size() );
//#else
//	spline_r_ = gsl_spline_alloc( gsl_interp_cspline, z_orig.size() );
//#endif
//	gsl_spline_init( spline_r_, z_spline, r_spline, z_orig.size() );
//	accel_i_ = gsl_interp_accel_alloc();
//#if GSL_MAJOR_VERSION > 1
//	spline_i_ = gsl_spline_alloc( gsl_interp_steffen, z_orig.size() );
//#else
//	spline_i_ = gsl_spline_alloc( gsl_interp_cspline, z_orig.size() );
//#endif
//	gsl_spline_init( spline_i_, z_spline, i_spline, z_orig.size() );
//
//	ierr = VecCreate( PETSC_COMM_SELF, psi );CHKERRQ(ierr);
//	ierr = VecSetSizes( *psi, PETSC_DECIDE, NZ_new );CHKERRQ(ierr);
//	ierr = VecSetFromOptions( *psi ); CHKERRQ(ierr);
//	ierr = VecSet( *psi, 0.0 );
//
//	for (ii = 0; ii < (PetscInt)NZ_new; ii++) {
//		tempval = gsl_spline_eval( spline_r_, z_new[ ii ], accel_r_ )
//				  + J * gsl_spline_eval( spline_i_, z_new[ ii ], accel_i_ );
//		ierr = VecSetValues( *psi, 1, &ii, &tempval, INSERT_VALUES );CHKERRQ(ierr);
//	}
//	ierr = VecAssemblyBegin( *psi );CHKERRQ(ierr);
//	ierr = VecAssemblyEnd( *psi );CHKERRQ(ierr);
//
//	gsl_spline_free( spline_r_ );
//	gsl_spline_free( spline_i_ );
//	gsl_interp_accel_free( accel_r_ );
//	gsl_interp_accel_free( accel_i_ );
//
//	delete [] z_spline;
//	delete [] r_spline;
//	delete [] i_spline;
//
//	return 1;
//}

int NCPA::EPadeSolver::build_operator_matrix_without_topography( 
	int NZvec, double *zvec, double k0, double h2, 
	std::complex<double> impedence_factor, std::complex<double> *n, size_t nqp, 
	int boundary_index, Mat *q ) {

	// Mat q;
	PetscInt Istart, Iend, col[3];
	PetscBool FirstBlock = PETSC_FALSE, LastBlock = PETSC_FALSE;
	PetscErrorCode ierr;
	PetscScalar value[3];
	PetscInt i;

	// Set up matrices
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZvec, NZvec, 3, NULL, q );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *q );CHKERRQ(ierr);

	// populate
	//double bnd_cnd = -1.0 / h2;    // @todo add hook for alternate boundary conditions
	//double bnd_cnd = -2.0 / h2;      // pressure release condition
	std::complex<double> bnd_cnd = (impedence_factor * std::sqrt( h2 ) - 1.0) / h2;
	double k02 = k0*k0;
	
	ierr = MatGetOwnershipRange(*q,&Istart,&Iend);CHKERRQ(ierr);
	if (Istart==0) FirstBlock=PETSC_TRUE;
    if (Iend==NZvec) LastBlock=PETSC_TRUE;
    value[0]=1.0 / h2 / k02; value[2]=1.0 / h2 / k02;
    for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
    		if (i < boundary_index)  {
    			value[ 0 ] = 0.0;
    			value[ 1 ] = 0.0;  
    			value[ 2 ] = 0.0;
    		} else if (i == boundary_index) {
    			value[ 0 ] = 0.0;
    			value[ 1 ] = bnd_cnd/k02 + (n[i]*n[i] - 1);
    			value[ 2 ] = 1.0 / h2 / k02;
    		} else {
    			value[ 0 ] = 1.0 / h2 / k02;
    			value[ 1 ] = -2.0 / h2 / k02 + (n[i]*n[i] - 1);
    			value[ 2 ] = 1.0 / h2 / k02;
    		}
		    col[0]=i-1; col[1]=i; col[2]=i+1;
		    ierr = MatSetValues(*q,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (LastBlock) {
		    i=NZvec-1; col[0]=NZvec-2; col[1]=NZvec-1;
		    value[ 0 ] = 1.0 / h2 / k02;
		    value[ 1 ] = -2.0/h2/k02 + (n[i]*n[i] - 1);
		    ierr = MatSetValues(*q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (FirstBlock) {
		    i=0; col[0]=0; col[1]=1; 
		    if (i < boundary_index)  {
    			value[ 0 ] = 0.0;
    			value[ 1 ] = 0.0;
    		} else {
    			value[ 0 ] = bnd_cnd/k02 + (n[i]*n[i] - 1);
    			value[ 1 ] = 1.0 / h2 / k02;
    		}
		    //value[0]=bnd_cnd/k02 + (n[i]*n[i] - 1); 
		    //value[1]=1.0 / h2 / k02;
		    ierr = MatSetValues(*q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(*q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    return 1;
}

double NCPA::EPadeSolver::check_ground_height_coincidence_with_grid( double *z, 
	size_t NZ, double tolerance, double z_ground ) {

	int closest_source_grid_point = (int)(NCPA::find_closest_index( z, NZ, z_ground ));
	if (std::fabs(z_ground - z[ closest_source_grid_point ]) < tolerance) {
		if (z_ground <= z[ closest_source_grid_point ]) {
			// ground is below grid point, move further down
			return z[ closest_source_grid_point ] - tolerance;
		} else {
			// ground is above grid point, move further up
			return z[ closest_source_grid_point ] + tolerance;
		}
	} else {
		return z_ground;
	}
}

int NCPA::EPadeSolver::solve_with_topography() {
	size_t i;
	std::complex<double> I( 0.0, 1.0 );
	PetscErrorCode ierr;
	PetscInt *indices;
	PetscScalar hank, *contents;
	Mat B, C, q;    // , q_starter;
	Mat *qpowers = PETSC_NULL, *qpowers_starter = PETSC_NULL;
	Vec psi_o, Bpsi_o; //, psi_temp;
	KSP ksp;
	// PC pc;

	// for turbulence, if needed
	double *mu_r, *mu_rpdr;
	std::vector<double> rand1, rand2;

	// set up z grid for flat ground.  When we add terrain we will need to move this inside
	// the range loop
	int profile_index;
	double minlimit, maxlimit;
	atm_profile_2d->get_maximum_altitude_limits( minlimit, maxlimit );
//	z_max = NCPA::min( z_max, minlimit );    // lowest valid top value
	z_max.set_value( NCPA::min( z_max.get(), minlimit ) );
	int ground_index = 0;
	std::complex<double> ground_impedence_factor( 0.0, 0.0 );
	if (r_max.get() > atm_profile_2d->get_maximum_valid_range() ) {
			atm_profile_2d->set_maximum_valid_range( r_max.get() );
		}

	// truncate multiprop file if needed
	if (write2d) {
		std::ofstream ofs( tag_filename(NCPAPROP_EPADE_PE_FILENAME_2D),
			std::ofstream::out | std::ofstream::trunc );
		ofs.close();
	}

//	// truncate 1-D if necessary
//	if (broadband) {
//		std::ofstream ofs( tag_filename(NCPAPROP_EPADE_PE_FILENAME_1D),
//			std::ofstream::out | std::ofstream::trunc );
//		ofs.close();
//	}

	// truncate topography if necessary
	if (write_topo) {
		std::ofstream ofs( tag_filename(NCPAPROP_EPADE_PE_FILENAME_TOPOGRAPHY),
			std::ofstream::out | std::ofstream::trunc );
		ofs.close();
	}

	double c0;
	for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm_profile_2d->first_profile();
		 it != atm_profile_2d->last_profile(); ++it) {
		if ( (*it)->contains_vector( "C0" ) ) {
			(*it)->convert_property_units( "C0", NCPAPROP_EPADE_PE_UNITS_C );
			(*it)->copy_vector_property( "C0", "_C0_" );
			c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
		} else {
			if ( (*it)->contains_vector("P") && (*it)->contains_vector("RHO") ) {
				(*it)->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO",
						NCPAPROP_EPADE_PE_UNITS_C );
				c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
			} else if ( (*it)->contains_vector("T") ) {
				(*it)->calculate_sound_speed_from_temperature( "_C0_", "T",
						NCPAPROP_EPADE_PE_UNITS_C );
				c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
			} else if ( (*it)->contains_vector( "CEFF" ) ) {
				c0 = atm_profile_2d->get( 0.0, "CEFF", z_ground );
			} else {
				throw std::runtime_error( "Cannot calculate static sound speed: None of CEFF, C0, T, or (P and RHO) are specified in atmosphere.");
			}
		}
	}

	// calculate/check z resolution
//	dz = param->getFloat( "dz_m" );
//	double f0 = f_vector[0];
	// @todo move this into f loop for future broadband use
	double dz;
	double lambda0 = c0 / f_vector[0];
	if (dz_requested.get() <= 0.0) {
		dz = lambda0 / 20.0;
		double nearestpow10 = std::pow( 10.0, (double)std::floor( (double)std::log10( dz ) ) );
		double factor = std::floor( dz / nearestpow10 );
		dz = nearestpow10 * factor;
		std::cout << "Setting dz to " << dz << " m" << std::endl;
	}
	if (dz > (c0 / f_vector[0] / 10.0) ) {
		std::ostringstream oss;
		oss << "Altitude resolution is too coarse.  Must be <= " << lambda0 / 10.0 << " meters.";
		throw std::runtime_error( oss.str() );
	}

	/* @todo move this into constructor as much as possible */
	z_bottom = -NCPAPROP_EPADE_PE_BASEMENT_THICKNESS;    // make this eventually depend on frequency
	z_bottom -= std::fmod( z_bottom, dz );
	// z_ground = atm_profile_2d->get( 0.0, "Z0" );
	z_ground = atm_profile_2d->get_interpolated_ground_elevation( 0.0 );
	NZ = ((int)std::floor((z_max.get() - z_bottom) / dz)) + 1;
	z = NCPA::zeros<double>( NZ );
	z_abs = NCPA::zeros<double>( NZ );
	indices = NCPA::zeros<PetscInt>( NZ );
	for (i = 0; i < NZ; i++) {
		z[ i ] = ((double)i) * dz + z_bottom;
		z_abs[ i ] = z[ i ];
		indices[ i ] = i;
	}
	double grid_tolerance =
		dz * NCPAPROP_EPADE_PE_GRID_COINCIDENCE_TOLERANCE_FACTOR;
	z_ground = check_ground_height_coincidence_with_grid(
		z, NZ, grid_tolerance, z_ground );
	double zs_abs = z_ground + z_source.get();
	// zs = NCPA::max( zs, z_ground );

	// define ground_index, which is J in @notes
	ground_index = (int)(NCPA::find_closest_index( z, NZ, z_ground ));
	if ( z[ ground_index ] < z_ground ) {
		ground_index++;
	}

	// adjust source height if it falls within 5% of a ground point
	int closest_source_grid_point =
		(int)(NCPA::find_closest_index( z, NZ, zs_abs ));
	if (std::fabs(zs_abs - z[ closest_source_grid_point ])
			< grid_tolerance) {
		zs_abs = z[ closest_source_grid_point ] + grid_tolerance;
		std::cout << "Adjusting source height to " << zs_abs-z_ground
			<< " m to avoid grid point singularity" << std::endl;
	}

	if (use_turbulence) {
		if (!random_turbulence) {
			// if (verbose) {
				std::cout << "Reading " << 2*turbulence_size
						  << " values from " << turbulence_file
						  << std::endl;
			// }
			std::ifstream rand_in( turbulence_file );
			if (!rand_in.good()) {
				throw std::runtime_error(
					"Error opening " + turbulence_file);
			}
			rand1.reserve( turbulence_size );
			for (i = 0; i < turbulence_size; i++) {
				rand_in >> rand1[ i ];
				if (!rand_in.good()) {
					std::ostringstream oss;
					oss << "Error reading turbulence numbers from "
						<< turbulence_file;
					throw std::runtime_error(oss.str());
				}
			}
			rand2.reserve( turbulence_size );
			for (i = 0; i < turbulence_size; i++) {
				rand_in >> rand2[ i ];
			}
			rand_in.close();
		}
	}
	
	// constants for now
	double h = dz;
	double h2 = h * h;
	double dr;

	// set up for source atmosphere
	double k0 = 0.0; c0 = 0.0;
	double *c = NCPA::zeros<double>( NZ );
	double *a_t = NCPA::zeros<double>( NZ );
	std::complex<double> *k = NCPA::zeros<std::complex<double>>( NZ );
	std::complex<double> *n = NCPA::zeros<std::complex<double>>( NZ );

	std::complex<double> *source = NCPA::zeros<std::complex<double>>( NZ );
	// @todo account for unspecified z_source being on the ground, which could be at a negative elevation instead of 0
	if (starter_type == NCPA::StarterType::SELF) {
		if (pointsource) {
			std::cout << "Generating point source at " << z_source.get() << "m" << std::endl;
			make_point_source( NZ, z, z_source.get() + z_ground, z_ground, source );
		} else {
			std::cout << "Reading line source from " << linesourcefile << std::endl;
			read_line_source_from_file( NZ, z, z_ground,
				linesourcefile, source );
		}

		// output source file for checking
		if (_write_source_function) {
			std::cout << "Writing source function to "
					<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_SOURCE) << std::endl;
			write_source(tag_filename(NCPAPROP_EPADE_PE_FILENAME_SOURCE), source, z, NZ);
		}
	}

	// write broadband header for testing
//	if (broadband) {
//		write_broadband_header( tag_filename(NCPAPROP_EPADE_PE_FILENAME_BROADBAND),
//			azi, NAz, f, NF, 1.0e8 );
//	}

	// freq and calc_az hold the current values of azimuth and frequency, respectively
	// these are used in the output routines, so make sure they get set correctly
	// whenever you change frequencies and azimuths
	for (std::vector<double>::const_iterator az = az_vector.cbegin(); az != az_vector.cend(); ++az) {
//	for (size_t azind = 0; azind < NAz; azind++) {
//		calc_az = azi[ azind ];
//		calc_az = *az;

		profile_index = -1;
		calculate_effective_sound_speed( atm_profile_2d, *az, "_CEFF_" );

		// atm_profile_2d->calculate_wind_component( "_WC_", "_WS_", "_WD_",
		// 	calc_az );
		// atm_profile_2d->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );

//		for (size_t freqind = 0; freqind < NF; freqind++) {
		for (std::vector<double>::const_iterator freq = f_vector.cbegin();
				freq != f_vector.cend(); ++freq) {

//			freq = f[ freqind ];
			std::cout << "Infrasound PE code at f = " << *freq << " Hz, azi = "
						<< *az << " deg" << std::endl;

			// calculate attenuation as a function of frequency if not externally supplied
//			if ( (!lossless) && (attnfile.length() == 0)) {
			if (attenuation_type == NCPA::AttenuationType::SUTHERLAND_BASS) {
				atm_profile_2d->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", *freq );
			}

			// Set up range vector
			if (nr_requested == 0) {
		  		dr = 340.0 / *freq;
				NR = (int)std::ceil( r_max.get() / dr );
		  	} else {
		  		NR = nr_requested;
		  		dr = r_max.get() / NR;
		  	}
		  	r = NCPA::zeros<double>( NR );
		  	zgi_r = NCPA::zeros<int>( NR );
		  	for (i = 0; i < NR; i++) {
		  		r[ i ] = ((double)(i+1)) * dr;
		  	}

		  	// Check ground elevation for coincidence with grid points
		  	z_ground = check_ground_height_coincidence_with_grid( z, NZ, 
		  		grid_tolerance,
		  		atm_profile_2d->get_interpolated_ground_elevation( 0.0 ) );

		  	// set up transmission loss matrix
			tl = NCPA::cmatrix( NZ, NR-1 );

			// calculate ground impedence (script A in notes in eq. 12)
			double rho0 = atm_profile_2d->get( 0.0, "RHO", z_ground );
			double lambBC = atm_profile_2d->get_first_derivative( 0.0, "RHO", z_ground ) / (2.0 * rho0);
			//lambBC = 0.0;
			if (user_ground_impedence_found) {
				ground_impedence_factor = *freq * I * 2.0 * PI * rho0 / user_ground_impedence + lambBC;
				std::cout << "Using user ground impedence of " << user_ground_impedence << std::endl;
			} else {
				ground_impedence_factor.real( lambBC );
				ground_impedence_factor.imag( 0.0 );
				std::cout << "Using default rigid ground with Lamb BC" << std::endl;
			}

			calculate_atmosphere_parameters( atm_profile_2d, NZ, z, 0.0, z_ground, attenuation_type,
				top_layer, *freq, use_topo, k0, c0, c, a_t, k, n );
			
			// calculate turbulence
			if (use_turbulence) {
				mu_r = NCPA::zeros<double>( NZ );
				mu_rpdr = NCPA::zeros<double>( NZ );
				setup_turbulence(rand1, rand2);
			}

			// build appropriate starter
			if (starter_type == NCPA::StarterType::SELF) {

				// for now build the non-topographic starter
				// revisit when time and funding permit
				size_t NZ_starter = NZ - ground_index;
				double *z_starter = NCPA::zeros<double>( NZ_starter ),
					   *c_starter = NCPA::zeros<double>( NZ_starter ),
					   *a_starter = NCPA::zeros<double>( NZ_starter );
				std::complex<double> *k_starter =
					NCPA::zeros<std::complex<double>>( NZ_starter );
				std::complex<double> *n_starter =
					NCPA::zeros<std::complex<double>>( NZ_starter );
				std::complex<double> *source_starter =
					NCPA::zeros<std::complex<double>>( NZ_starter );

				double k0_starter, c0_starter;
				size_t ii;
				for (ii = 0; ii < NZ_starter; ii++) {
					z_starter[ ii ] = z[ ii + ground_index ];
					source_starter[ ii ] = source[ ii + ground_index ];
				}

				calculate_atmosphere_parameters( atm_profile_2d, NZ_starter, z_starter,
					0.0, z_ground, attenuation_type, top_layer, *freq, false, k0_starter,
					c0_starter, c_starter, a_starter, k_starter, n_starter );

				Mat q_starter = PETSC_NULL;
				Mat *qpowers_starter = PETSC_NULL;
				build_operator_matrix_without_topography( NZ_starter, z_starter, 
					k0_starter, h2, ground_impedence_factor, n_starter, pade_order+1, 0,
					&q_starter );
				create_matrix_polynomial( pade_order+1, &q_starter, &qpowers_starter );
				get_starter_self( NZ_starter, z_starter, source_starter, k0_starter,
					qpowers_starter, pade_order, &psi_o );
//				outputVec( psi_o, z_starter, NZ_starter,
//						"starter.prespline.dat" );

				// now interpolate calculated starter to actual Z vector
//				std::deque< double > z_spline, r_spline, i_spline;
				PetscScalar *psi_orig = NCPA::zeros<PetscScalar>( NZ_starter );
				PetscScalar *psi_new = NCPA::zeros<PetscScalar>( NZ );
				PetscInt *starter_indices = NCPA::index_vector<PetscInt>( NZ_starter );

				ierr = VecGetValues( psi_o, NZ_starter, starter_indices, psi_orig );CHKERRQ(ierr);
//				for (ii = 0; ii < NZ_starter; ii++) {
//					z_spline.push_back( z_starter[ ii ] );
//					r_spline.push_back( psi_orig[ ii ].real() );
//					i_spline.push_back( psi_orig[ ii ].imag() );
//				}
				ierr = VecDestroy( &psi_o );
				interpolate_complex( NZ_starter, z_starter, psi_orig, NZ, z, psi_new );
				ierr = VecCreate( PETSC_COMM_SELF, &psi_o );CHKERRQ(ierr);
				ierr = VecSetSizes( psi_o, PETSC_DECIDE, NZ );CHKERRQ(ierr);
				ierr = VecSetFromOptions( psi_o ); CHKERRQ(ierr);
				ierr = VecSet( psi_o, 0.0 );
				PetscInt *z_indices = NCPA::index_vector<PetscInt>( NZ );
				ierr = VecSetValues( psi_o, NZ, z_indices, psi_new, INSERT_VALUES );CHKERRQ(ierr);
				ierr = VecAssemblyBegin( psi_o );CHKERRQ(ierr);
				ierr = VecAssemblyEnd( psi_o );CHKERRQ(ierr);

				// clean up temp variables
				delete [] z_starter;
				delete [] k_starter;
				delete [] n_starter;
				delete [] c_starter;
				delete [] a_starter;
				delete [] starter_indices;
				delete [] z_indices;
				delete [] psi_orig;
				delete [] psi_new;
				ierr = MatDestroy( &q_starter );CHKERRQ(ierr);
				delete_matrix_polynomial( pade_order+1, &qpowers_starter );


				// set up for future calculations
				build_operator_matrix_with_topography( atm_profile_2d, NZ, z, 0.0, k, k0, 
					h2, z_ground, ground_impedence_factor, n, ground_index, PETSC_NULL, 
					&q, true );
				create_matrix_polynomial( pade_order+1, &q, &qpowers );
				
				ierr = MatDestroy( &q );CHKERRQ(ierr);
			} else if (starter_type == NCPA::StarterType::GAUSSIAN) {
				// qpowers_starter = qpowers;
				build_operator_matrix_with_topography( atm_profile_2d, NZ, z, 0.0, k, k0, 
					h2, z_ground, ground_impedence_factor, n, ground_index, PETSC_NULL, 
					&q, true );
				create_matrix_polynomial( pade_order+1, &q, &qpowers );
				ierr = MatDestroy( &q );CHKERRQ(ierr);
				get_starter_gaussian( NZ, z, z_source.get()+z_ground, k0,
					ground_index, &psi_o );
			} else if (starter_type == NCPA::StarterType::USER) {
				build_operator_matrix_with_topography( atm_profile_2d, NZ, z, 0.0, k, k0, 
					h2, z_ground, ground_impedence_factor, n, ground_index, PETSC_NULL, 
					&q, true );
				create_matrix_polynomial( pade_order+1, &q, &qpowers );
				ierr = MatDestroy( &q );CHKERRQ(ierr);
				get_starter_user( starter_filename_in, NZ, z, &psi_o );
			} else {
				std::cerr << "Unrecognized starter type:" << std::endl;
				exit(0);
			}

			if (write_starter) {
				std::cout << "Outputting starter..." << std::endl;
				NCPA::outputVec( psi_o, z, NZ, tag_filename(NCPAPROP_EPADE_PE_FILENAME_STARTER) );
			}

			std::cout << "Finding ePade coefficients..." << std::endl;
			std::vector< PetscScalar > P, Q;
			std::vector< PetscScalar > taylor = taylor_exp_id_sqrt_1pQ_m1( 2*pade_order, k0*dr );
			calculate_pade_coefficients( &taylor, pade_order, pade_order+1, &P, &Q );
			generate_polymatrices( qpowers, pade_order, NZ, P, Q, &B, &C );

			std::cout << "Marching out field..." << std::endl;
			ierr = VecDuplicate( psi_o, &Bpsi_o );CHKERRQ(ierr);
			contents = NCPA::zeros<PetscScalar>( NZ );

			ierr = KSPCreate( PETSC_COMM_SELF, &ksp );CHKERRQ(ierr);
			ierr = KSPSetOperators( ksp, C, C );CHKERRQ(ierr);
			ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
			for (size_t ir = 0; ir < (NR-1); ir++) {

				double rr = r[ ir ];
				z_ground = check_ground_height_coincidence_with_grid( z, NZ, 
		  			grid_tolerance,
		  			atm_profile_2d->get_interpolated_ground_elevation( rr ) );
				calculate_atmosphere_parameters( atm_profile_2d, NZ, z, rr, z_ground,
					attenuation_type, top_layer, *freq, use_topo, k0, c0, c, a_t, k, n );
				Mat last_q;
				ierr = MatConvert( qpowers[0], MATSAME, MAT_INITIAL_MATRIX, &last_q );CHKERRQ(ierr);
				delete_matrix_polynomial( pade_order+1, &qpowers );
				
				build_operator_matrix_with_topography( atm_profile_2d, NZ, z, rr, k, 
					k0, h2, z_ground, ground_impedence_factor, n, ground_index, 
					last_q, &q );
				create_matrix_polynomial( pade_order+1, &q, &qpowers );
				
				ierr = MatDestroy( &last_q );CHKERRQ(ierr);
				ierr = MatDestroy( &q );CHKERRQ(ierr);

				taylor = taylor_exp_id_sqrt_1pQ_m1( 2*pade_order, k0*dr );
				calculate_pade_coefficients( &taylor, pade_order, pade_order+1, &P, &Q );
				ierr = MatZeroEntries( B );CHKERRQ(ierr);
				ierr = MatZeroEntries( C );CHKERRQ(ierr);
				generate_polymatrices( qpowers, pade_order, NZ, P, Q, &B, &C );
				
				// apply turbulence
				ierr = VecGetValues( psi_o, NZ, indices, contents );
				if (use_turbulence) {
					// update ground index
					ground_index = (int)(NCPA::find_closest_index( z, NZ, z_ground ));
					if ( z[ ground_index ] < z_ground ) {
						ground_index++;
					}
					if (ir == 0) {
						// calculate first step
						calculate_turbulence( rr, NZ, z, k0, ground_index,
							mu_r );
					} else {
						std::memcpy( mu_r, mu_rpdr, NZ*sizeof(double) );
					}
					calculate_turbulence( rr + dr, NZ, z, k0, ground_index,
						mu_rpdr );

					// apply the turbulent fluctuations.  Do this inside
					// the if() because we need to keep these modifications
					// to psi_o, as opposed to the scaling by the Hankel
					// function below
					for (i = 0; i < NZ; i++) {
						contents[ i ] *= std::exp( I * k0 * dr * 0.5 *
							(mu_r[ i ] + mu_rpdr[ i ]) );
					}

					// store the modified field
					ierr = VecSetValues( psi_o, NZ, indices, contents,
						INSERT_VALUES );CHKERRQ(ierr);
					ierr = VecAssemblyBegin( psi_o );CHKERRQ(ierr);
					ierr = VecAssemblyEnd( psi_o );CHKERRQ(ierr);
				}

				hank = std::sqrt( 2.0 / ( PI * k0 * rr ) ) * exp( I * ( k0 * rr - PI/4.0 ) );
				for (i = 0; i < NZ; i++) {
					tl[ i ][ ir ] = contents[ i ] * hank;
				}

				// make sure the receiver height is above ground
				double z0g = z_ground + z_receiver.get();
				// z0g = NCPA::max( z0g, zr );
				zgi_r[ ir ] = (int)(NCPA::find_closest_index( z, NZ, z0g ));
				while ( z[ zgi_r[ ir ] ] <= z_ground ) {
//					std::cout << "Adjusting reported z for r[" << ir << "] = "
//							<< r[ ir ] << " from "
//							<< z[ zgi_r[ ir ] ] << " to " << z[ zgi_r[ ir ]+1 ]
//							<< " for z_g = " << z_ground << std::endl;
					zgi_r[ ir ]++;
				}
//				zgi_r[ir]++;
				
				if ( std::fmod( rr, 1.0e5 ) < dr) {
					std::cout << " -> Range " << rr/1000.0 << " km" << std::endl;
				}

				// Set up the system
				// C * psi_o_next = B * psi_o
				ierr = MatMult( B, psi_o, Bpsi_o );CHKERRQ(ierr);
				ierr = KSPSetOperators( ksp, C, C );CHKERRQ(ierr);
				ierr = KSPSolve( ksp, Bpsi_o, psi_o );CHKERRQ(ierr);
			}
			std::cout << "Stopped at range " << r[ NR-1 ]/1000.0 << " km" << std::endl;

			

			if (az_vector.size() > 1) {
				if (write1d) {
					std::cout << "Writing 1-D output to "
						<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP) << std::endl;
					output1DTL( tag_filename(NCPAPROP_EPADE_PE_FILENAME_MULTIPROP), *az, true );
				}
			} else { 
				if (write1d) {
					std::cout << "Writing 1-D output to "
						<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_1D) << std::endl;
					output1DTL( tag_filename(NCPAPROP_EPADE_PE_FILENAME_1D), *az, false );
				}
				if (write2d) {
					std::cout << "Writing 2-D output to "
						<< tag_filename(NCPAPROP_EPADE_PE_FILENAME_2D) << std::endl;
					output2DTL( tag_filename(NCPAPROP_EPADE_PE_FILENAME_2D) );
				}
			}

			// write broadband body for testing
//			if (broadband) {
//				write_broadband_results( tag_filename(NCPAPROP_EPADE_PE_FILENAME_BROADBAND),
//					calc_az, freq, r, NR, z_abs, NZ, tl, 1.0e8 );
//			}

			if (write_atmosphere) {
				std::cout << "Writing source atmosphere to "
						<< tag_filename("atm_profile.pe") << std::endl;
				std::vector<std::string> keylist;
				keylist.push_back( "U" );
				keylist.push_back( "V" );
				keylist.push_back( "T" );
				keylist.push_back( "RHO" );
				keylist.push_back( "P" );
				keylist.push_back( "_C0_" );
				keylist.push_back( "_CEFF_" );
				std::ofstream atmout( tag_filename( "atm_profile.pe" ) );
				atm_profile_2d->print_atmosphere( keylist, 0.0, "Z", atmout );
				atmout.close();
			}

			if (use_turbulence) {
				delete [] mu_r;
				delete [] mu_rpdr;
				cleanup_turbulence();
			}
			
			std::cout << std::endl;

			delete_matrix_polynomial( pade_order+1, &qpowers );
			if (starter_type == NCPA::StarterType::SELF) {
				delete_matrix_polynomial( pade_order+1, &qpowers_starter );
			}
//			if (attnfile.length() == 0) {
			if (attenuation_type != NCPA::AttenuationType::USER) {
				atm_profile_2d->remove_property( "_ALPHA_" );
			}
			delete [] r;
			delete [] zgi_r;
			NCPA::free_cmatrix( tl, NZ, NR-1 );
		}

		if (write_topo) {
			write_topography(
				tag_filename(NCPAPROP_EPADE_PE_FILENAME_TOPOGRAPHY),
				*az, r_max.get(), 1000.0 );
		}
		
		atm_profile_2d->remove_property( "_CEFF_" );
		// atm_profile_2d->remove_property( "_WC_" );
	}

	ierr = MatDestroy( &B );       CHKERRQ(ierr);
	ierr = MatDestroy( &C );       CHKERRQ(ierr);
	ierr = VecDestroy( &psi_o );   CHKERRQ(ierr);
	ierr = VecDestroy( &Bpsi_o );  CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp );     CHKERRQ(ierr);
	
	delete [] k;
	delete [] n;
	delete [] c;
	delete [] a_t;
	delete [] contents;
	delete [] indices;
	delete [] z;
	delete [] z_abs;

	return 1;
}

// create a vector of powers of a given matrix
int NCPA::EPadeSolver::create_matrix_polynomial( size_t nterms, const Mat *Q, Mat **qpowers ) {

	PetscErrorCode ierr;
	PetscInt i;

	if ((*qpowers) != PETSC_NULL) {
		delete_matrix_polynomial( nterms, qpowers );
	}

	*qpowers = new Mat[ nterms ];
	ierr = MatConvert( *Q, MATSAME, MAT_INITIAL_MATRIX, *qpowers );CHKERRQ(ierr);
	for (i = 1; i < (PetscInt)nterms; i++) {
		ierr = MatCreate( PETSC_COMM_SELF, (*qpowers) + i );CHKERRQ(ierr);
		ierr = MatSetFromOptions( (*qpowers)[ i ] );CHKERRQ(ierr);
		ierr = MatMatMult( (*qpowers)[i-1], (*qpowers)[0], MAT_INITIAL_MATRIX, PETSC_DEFAULT, 
			(*qpowers) + i );CHKERRQ(ierr);
	}

	return 1;
}

// clean up a vector of powers of a matrix
int NCPA::EPadeSolver::delete_matrix_polynomial( size_t nterms, Mat **qpowers ) {
	PetscErrorCode ierr;
	if ((*qpowers) != PETSC_NULL) {
		for (size_t i = 0; i < nterms; i++) {
			ierr = MatDestroy( (*qpowers) + i ); CHKERRQ(ierr);
		}
		delete [] *qpowers;
		*qpowers = PETSC_NULL;
	}

	return 1;
}

// Calculate and return k0, c0, c, a, k, and n
void NCPA::EPadeSolver::calculate_atmosphere_parameters(
	NCPA::Atmosphere2D *atm, int NZvec, double *z_vec,
	double r, double z_g, NCPA::AttenuationType attntype, bool use_top_layer, double freq, bool absolute,
	double &k0, double &c0, double *c_vec, double *a_vec, std::complex<double> *k_vec, 
	std::complex<double> *n_vec ) {

	std::complex<double> I( 0.0, 1.0 );

	std::memset( c_vec, 0, NZvec * sizeof(double) );
	std::memset( a_vec, 0, NZvec * sizeof(double) );
//	std::memset( k_vec, 0, NZvec * sizeof( std::complex< double > ) );
//	std::memset( n_vec, 0, NZvec * sizeof( std::complex< double > ) );
	std::fill( k_vec, k_vec+NZvec, std::complex<double>{} );
	std::fill( n_vec, n_vec+NZvec, std::complex<double>{} );

	// z_vec is relative to ground
	if (absolute) {
		fill_atm_vector_absolute( atm, r, NZvec, z_vec, "_CEFF_", c_underground, c_vec );
	} else {
		fill_atm_vector_relative( atm, r, NZvec, z_vec, "_CEFF_", z_g, c_vec );
	}
	c0 = atm->get( r, "_CEFF_", z_g );

//	if (!use_lossless) {
	if (attntype != NCPA::AttenuationType::NONE) {
		if (absolute) {
			fill_atm_vector_absolute( atm, r, NZvec, z_vec, "_ALPHA_", 0.0, a_vec );
		} else {
			fill_atm_vector_relative( atm, r, NZvec, z_vec, "_ALPHA_", z_g, a_vec );
		}
	}
	double *abslayer = NCPA::zeros<double>( NZvec );
	if (use_top_layer) {
		double tlt = top_layer_thickness_m;
		if (tlt < 0.0) {
			tlt = NCPA::min<double>( c0 / freq, 5000.0 );
		}
		absorption_layer( tlt, z_vec, NZvec, abslayer );
		// std::ofstream absout("abslayer_2d.dat");
		// NCPA::print_2_columns<double,double>( absout, NZ, z, abslayer );
		// absout.close();
	}

	k0 = 2.0 * PI * freq / c0;
	
	// Set up vectors
	for (int i = 0; i < NZvec; i++) {
		double rho, drho, ddrho;
		if (absolute && (z_vec[i] < z_g)) {
			k_vec[ i ] = 0.0;    // k == 0 below the ground
		} else {
			rho = atm->get( r, "RHO", z_vec[ i ] );
			drho = atm->get_first_derivative( r, "RHO", z_vec[ i ] );
			ddrho = atm->get_second_derivative( r, "RHO", z_vec[ i ] );
			k_vec[ i ] = std::sqrt(
							std::pow( 2.0 * PI * freq / c_vec[ i ], 2.0 )
							- 0.75 * std::pow( drho / rho, 2.0 )
							+ 0.5 * ddrho / rho
						) + (a_vec[ i ] + abslayer[ i ]) * I;
			//k_vec[ i ] = 2.0 * PI * freq / c_vec[ i ] + a_vec[ i ] * I;
		}
		n_vec[ i ] = k_vec[ i ] / k0;
	}
}

void NCPA::EPadeSolver::fill_atm_vector_relative( NCPA::Atmosphere2D *atm, double range, int NZvec, double *zvec, 
	std::string key, double groundheight, double *vec ) {

	for (int i = 0; i < NZvec; i++) {
		vec[i] = atm->get( range, key, zvec[i] + groundheight );
	}
}

void NCPA::EPadeSolver::fill_atm_vector_absolute( NCPA::Atmosphere2D *atm, double range, int NZvec, double *zvec, 
	std::string key, double fill_value, double *vec ) {

	double zmin = atm->get_interpolated_ground_elevation( range ); 

	for (int i = 0; i < NZvec; i++) {
		if (zvec[i] < zmin) {
			vec[i] = fill_value;
		} else {
			vec[i] = atm->get( range, key, zvec[i] );
		}
	}
}


int NCPA::EPadeSolver::generate_polymatrix( Mat *qpowers, size_t qpowers_size, size_t NZ,
			std::vector< std::complex< double > > &T, Mat *B ) {

	PetscErrorCode ierr;
	PetscInt Istart, Iend, i;
	PetscScalar value;

	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*T.size()+1, NULL, B );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *B );CHKERRQ(ierr);

	// start B off as T[0]
	ierr = MatGetOwnershipRange(*B,&Istart,&Iend);CHKERRQ(ierr);
	value = T[0];
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( *B, 1, &i, 1, &i, &value, INSERT_VALUES );CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    PetscInt nterms = NCPA::min( T.size(), (size_t)qpowers_size );
	for (i = 1; i < (PetscInt)nterms; i++) {
		ierr = MatAXPY( *B, T[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	return 1;
}

int NCPA::EPadeSolver::sum_scaled_matrix_polynomial_terms( Mat *qpowers, int qpowers_size,
	int NZ, std::vector< std::complex< double > > &T, Mat *B ) {

	PetscErrorCode ierr;
	PetscInt Istart, Iend, i;
	PetscScalar value;

	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*T.size()+1, NULL, B );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *B );CHKERRQ(ierr);

	// start B off as T[0]
	ierr = MatGetOwnershipRange(*B,&Istart,&Iend);CHKERRQ(ierr);
	value = T[0];
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( *B, 1, &i, 1, &i, &value, INSERT_VALUES );CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    PetscInt nterms = NCPA::min( T.size(), (size_t)qpowers_size );
	for (i = 1; i < (PetscInt)nterms; i++) {
		// std::cout << "Multiplying M ^ " << i << " by " << T[ i ] << std::endl;
		ierr = MatAXPY( *B, T[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	return 1;
}


int NCPA::EPadeSolver::generate_polymatrices( Mat *qpowers, size_t npade, size_t NZ,
	std::vector< std::complex< double > > &P, std::vector< std::complex< double > > &Q,
	Mat *B, Mat *C ) {

	PetscErrorCode ierr;
	PetscInt Istart, Iend, i;
	PetscScalar value;

	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade-1, NULL, B );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *B );CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade+1, NULL, C );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *C );CHKERRQ(ierr);

	ierr = MatGetOwnershipRange(*B,&Istart,&Iend);CHKERRQ(ierr);
	// by definition Q[0] is 1.  It so happens that P[0] is also 1, but this is not 
	// guaranteed.
	// @todo generalize this fot the case where P[0] != 1
	value = 1.0;
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( *B, 1, &i, 1, &i, &value, INSERT_VALUES );CHKERRQ(ierr);
	}
	ierr = MatGetOwnershipRange( *C, &Istart, &Iend );CHKERRQ(ierr);
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( *C, 1, &i, 1, &i, &value, INSERT_VALUES );CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(*C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    // std::cout << "Denominator:" << std::endl;
	for (i = 1; i < (PetscInt)(Q.size()); i++) {
		// std::cout << "Multiplying Q ^ " << i << " ] by " << Q[ i ] << std::endl;
		ierr = MatAXPY( *C, Q[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	// std::cout << "Numerator:" << std::endl;
	for (i = 1; i < (PetscInt)(P.size()); i++) {
		// std::cout << "Multiplying P ^ " << i << " by " << P[ i ] << std::endl;
		ierr = MatAXPY( *B, P[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	return 1;
}

int NCPA::EPadeSolver::build_operator_matrix_with_topography( NCPA::Atmosphere2D *atm, 
	int NZvec, double *zvec, double r, std::complex<double> *k, double k0, double h2, 
	double z_s, std::complex<double> impedence_factor, std::complex<double> *n, 
	int boundary_index, const Mat &last_q, Mat *q, bool starter ) {

	//Mat q;
	PetscInt Istart, Iend, *col, *indices;
	PetscInt *nonzeros;
	PetscBool FirstBlock = PETSC_FALSE, LastBlock = PETSC_FALSE;
	PetscErrorCode ierr;
	PetscScalar value[3], *rowDiff;
	PetscInt i, j;

	// Calculate parameters
	double h = std::sqrt( h2 );
	double J_s = (z_s - zvec[0]) / h;
	int Ji = (int)(NCPA::find_closest_index( zvec, NZvec, z_s ));  // first index above ground
	if ( zvec[ Ji ] < z_s ) {
		Ji++;
	}
	double dJ = (double)Ji;
	//double dJ = zvec[ Ji ] / h;

	// number of nonzero values
	nonzeros = NCPA::single_valued_vector<PetscInt>( NZvec, 3 );
	indices  = NCPA::index_vector<PetscInt>( NZvec );
	col      = NCPA::zeros<PetscInt>( NZvec );
	nonzeros[ 0 ] = 2;
	nonzeros[ NZvec-1 ] = 2;


	// calculate intermediate variables as shown in notes
	double rho_a, rho_b, Gamma;
	double Anom, Bnom, denom;
	double s_A, s_B, a, b, c, alpha, beta, gamma;

	rho_a = atm->get( r, "RHO", z_s );
	Gamma = 0.5 * atm->get_first_derivative( r, "RHO", z_s ) / rho_a;
	rho_b = RHO_B;
	Anom = (1.0 / rho_a) * (1.0 / (dJ - J_s));
	if (starter) {
		Bnom = 0;    // for starter, rho_b = inf
	} else {
		Bnom = (1.0 / rho_b) * (1.0 / (J_s - dJ + 1.0));
	}
	denom = Anom + Bnom - (h * Gamma);
	s_A = Anom / denom;
	s_B = Bnom / denom;
	a   = s_B / (dJ - J_s);
	b   = (s_A - 2.0) / (dJ - J_s);
	c   = 1.0 / (dJ - J_s);
	alpha = 1.0 / (J_s - dJ + 1.0);
	beta  = (s_B - 2.0) / (J_s - dJ + 1.0);
	gamma = s_A / (J_s - dJ + 1.0);

	// Calculate matrix ratio representation of sqrt(1+Q)
	rowDiff = NCPA::zeros<PetscScalar>( NZvec );
	if (last_q != PETSC_NULL) {
		PetscScalar I( 0.0, 1.0 ), *rowAbove, *rowBelow;
		PetscScalar M = I * k0 * atm->get_interpolated_ground_elevation_first_derivative( r ) * h / denom;
		Vec vecAbove, vecBelow;
		PetscInt num_nonzeros;
		approximate_sqrt_1pQ( NZvec, &last_q, Ji, &vecBelow, &vecAbove, &num_nonzeros );
		ierr = VecScale( vecBelow, M );CHKERRQ(ierr);
		ierr = VecScale( vecAbove, M );CHKERRQ(ierr);
		
		// get the Ji'th and (Ji-1)'th rows of the M matrix
		rowBelow = NCPA::zeros<PetscScalar>( NZvec );
		ierr = VecGetValues( vecBelow, NZvec, indices, rowBelow );CHKERRQ(ierr);
		nonzeros[ Ji-1 ] = num_nonzeros;

		rowAbove = NCPA::zeros<PetscScalar>( NZvec );
		ierr = VecGetValues( vecAbove, NZvec, indices, rowAbove );CHKERRQ(ierr);
		nonzeros[ Ji ] = num_nonzeros;
		
		for (i = 0; i < NZvec; i++) {
			rowDiff[ i ] = rowAbove[ i ] - rowBelow[ i ];
		}

		delete [] rowAbove;
		delete [] rowBelow;
		ierr = VecDestroy( &vecAbove );CHKERRQ(ierr);
		ierr = VecDestroy( &vecBelow );CHKERRQ(ierr);
	} 

	// Set up matrices
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZvec, NZvec, 0, nonzeros, q );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *q );CHKERRQ(ierr);
	

	// populate
	double k02 = k0*k0;
	
	// If this process is being split over processors, we need to check to see
	// if this particular instance contains the first or last rows, because those
	// get filled differently
	ierr = MatGetOwnershipRange(*q,&Istart,&Iend);CHKERRQ(ierr);

	// Does this instance contain the first row?
	if (Istart == 0) {
		FirstBlock=PETSC_TRUE;
	}

	// Does this instance contain the last row?
    if (Iend==(PetscInt)NZ) {
    	LastBlock=PETSC_TRUE;
    }

    // iterate over block.  If this instance contains the first row, leave that one
    // for later, same for if this instance contains the last row.
    PetscScalar *Drow = NCPA::zeros<PetscScalar>( NZvec );
    for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {

		// set column numbers.  Since the matrix Q is tridiagonal (because input 
		// matrix D is tridiagonal and K is diagonal), column indices are
		// i-1, i, i+1
    	col[ 0 ] = i-1;
    	col[ 1 ] = i;
    	col[ 2 ] = i+1;

    	// Set values.  This will be the same unless we're at the indices immediately
    	// below or above the ground surface
    	if (i == (Ji-1)) {

    		if (last_q != PETSC_NULL) {
	    		// this is the alpha, beta, gamma row
	    		std::memcpy( col, indices, NZvec*sizeof(PetscInt) );
	    		for (j = 0; j < NZvec; j++) {
	    			Drow[ j ] = -rowDiff[ j ] / ( h2 * (J_s - dJ + 1.0));
	    		}

	    		Drow[ Ji-2 ] += alpha;
	    		Drow[ Ji-1 ] += beta;
	    		Drow[  Ji  ] += gamma;
	    	} else {   // M == 0
	    		Drow[0] = alpha;
	    		Drow[1] = beta;
	    		Drow[2] = gamma;
	    	}
    		
    	} else if (i == Ji) {
    		// this is the a, b, c row
    		if (last_q != PETSC_NULL) {
	    		std::memcpy( col, indices, NZvec*sizeof(PetscInt) );
	    		for (j = 0; j < NZvec; j++) {
	    			Drow[ j ] = -rowDiff[ j ] / ( h2 * (dJ - J_s));
	    		}

	    		Drow[ Ji-1 ] += a;
	    		Drow[  Ji  ] += b;
	    		Drow[ Ji+1 ] += c;
	    	} else {
	    		Drow[0] = a;
	    		Drow[1] = b;
	    		Drow[2] = c;
	    	}
    	} else {
    		Drow[0] = 1.0;
    		Drow[1] = -2.0;
    		Drow[2] = 1.0;
    	}

    	for (j = 0; j < nonzeros[ i ]; j++) {
    		if (col[j] == i) {
    			Drow[ j ] = ( (Drow[ j ] / h2) + k[ i ]*k[ i ] - k02 ) / k02;
    		} else {
    			Drow[ j ] /= (h2 * k02);
    		}
    	}
    	ierr = MatSetValues(*q,1,&i,nonzeros[ i ],col,Drow,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (LastBlock) {
		    i=NZ-1; col[0]=NZ-2; col[1]=NZ-1;
		    value[ 0 ] = 1.0 / h2 / k02;
		    //value[ 1 ] = -2.0/h2/k02 + (n[i]*n[i] - 1);
		    value[ 1 ] = ( (-2.0 / h2) + k[ i ]*k[ i ] - k02 ) / k02;
		    ierr = MatSetValues(*q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (FirstBlock) {
		    i=0; col[0]=0; col[1]=1; 
		    if (i == (Ji-1))  {
		    	if (last_q != PETSC_NULL) {
			    	std::memcpy( col, indices, NZvec*sizeof(PetscInt) );
			    	for (j = 0; j < NZvec; j++) {
		    			Drow[ j ] = -rowDiff[ j ] / ( h2 * (J_s - dJ + 1.0));
		    		}
		    		// std::memcpy( Drow, rowBelow, NZvec*sizeof(PetscScalar) );

		    		Drow[ 0 ] += beta;
		    		Drow[ 1 ] += gamma;
		    	} else {
		    		Drow[ 0 ] = beta;
		    		Drow[ 1 ] = gamma;
		    	}
    		} else if (i == Ji) {
    			if (last_q != PETSC_NULL) {
	    			std::memcpy( col, indices, NZvec*sizeof(PetscInt) );
	    			for (j = 0; j < NZvec; j++) {
		    			Drow[ j ] = -rowDiff[ j ] / ( h2 * (dJ - J_s));
		    		}
		    		// std::memcpy( Drow, rowBelow, NZvec*sizeof(PetscScalar) );

		    		Drow[ 0 ] += b;
		    		Drow[ 1 ] += c;
		    	} else {
		    		Drow[ 0 ] = b;
		    		Drow[ 1 ] = c;
		    	}
    		} else {
    			Drow[0] = -2.0;
	    		Drow[1] = 1.0;
    		}
    		for (j = 0; j < nonzeros[ i ]; j++) {
	    		if (col[j] == i) {
	    			Drow[ j ] = ( (Drow[ j ] / h2) + k[ i ]*k[ i ] - k02 ) / k02;
	    		} else {
	    			Drow[ j ] /= (h2 * k02);
	    		}
	    	}

    		ierr = MatSetValues(*q,1,&i,nonzeros[i],col,Drow,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(*q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    delete [] nonzeros;
    delete [] indices;
    delete [] col;
    delete [] Drow;
    delete [] rowDiff;
	return 1;
}

int NCPA::EPadeSolver::approximate_sqrt_1pQ( int NZvec, const Mat *Q, PetscInt Ji, Vec *vecBelow, Vec *vecAbove, PetscInt *nonzeros ) {

	PetscErrorCode ierr;
	

	/* for order (0,0) */

	// const PetscScalar diag = 1;
	// const PetscInt diagIndBelow = Ji-1;
	// const PetscInt diagIndAbove = Ji;
	// ierr = VecCreate( PETSC_COMM_SELF, vecBelow );CHKERRQ(ierr);
	// ierr = VecSetType( *vecBelow, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( *vecBelow );CHKERRQ(ierr);
	// ierr = VecSetSizes( *vecBelow, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecSet( *vecBelow, 0 );CHKERRQ(ierr);
	// ierr = VecSetValues( *vecBelow, 1, &diagIndBelow, &diag, INSERT_VALUES );CHKERRQ(ierr);
	// ierr = VecCreate( PETSC_COMM_SELF, vecAbove );CHKERRQ(ierr);
	// ierr = VecSetType( *vecAbove, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( *vecAbove );CHKERRQ(ierr);
	// ierr = VecSetSizes( *vecAbove, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecSet( *vecAbove, 0 );CHKERRQ(ierr);
	// ierr = VecSetValues( *vecBelow, 1, &diagIndAbove, &diag, INSERT_VALUES );CHKERRQ(ierr);
	// ierr = VecAssemblyBegin( *vecBelow );CHKERRQ(ierr);
	// ierr = VecAssemblyEnd( *vecBelow );CHKERRQ(ierr);
	// ierr = VecAssemblyBegin( *vecAbove );CHKERRQ(ierr);
	// ierr = VecAssemblyEnd( *vecBelow );CHKERRQ(ierr);
	// *nonzeros = 3;
	// return 1;

	/* For order (1,0) */
	
	PetscInt nvals;
	const PetscInt *indices;
	Mat halfQ;
	const PetscScalar *values;
	const PetscScalar diag = 1;
	const PetscInt diagIndBelow = Ji-1;
	const PetscInt diagIndAbove = Ji;
	
	ierr = VecCreate( PETSC_COMM_SELF, vecBelow );CHKERRQ(ierr);
	ierr = VecSetType( *vecBelow, VECSEQ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( *vecBelow );CHKERRQ(ierr);
	ierr = VecSetSizes( *vecBelow, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	ierr = VecSet( *vecBelow, 0 );CHKERRQ(ierr);
	ierr = VecCreate( PETSC_COMM_SELF, vecAbove );CHKERRQ(ierr);
	ierr = VecSetType( *vecAbove, VECSEQ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( *vecAbove );CHKERRQ(ierr);
	ierr = VecSetSizes( *vecAbove, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	ierr = VecSet( *vecAbove, 0 );CHKERRQ(ierr);

	ierr = MatDuplicate( *Q, MAT_COPY_VALUES, &halfQ );CHKERRQ(ierr);
	ierr = MatScale( halfQ, 0.5 );
	ierr = MatGetRow( halfQ, Ji-1, &nvals, &indices, &values );CHKERRQ(ierr);
	ierr = VecSetValues( *vecBelow, nvals, indices, values, INSERT_VALUES );CHKERRQ(ierr);
	ierr = VecSetValues( *vecBelow, 1, &diagIndBelow, &diag, ADD_VALUES );CHKERRQ(ierr);
	*nonzeros = nvals;
	ierr = MatRestoreRow( halfQ, Ji-1, &nvals, &indices, &values );CHKERRQ(ierr);

	ierr = MatGetRow( halfQ, Ji, &nvals, &indices, &values );CHKERRQ(ierr);
	ierr = VecSetValues( *vecAbove, nvals, indices, values, INSERT_VALUES );CHKERRQ(ierr);
	ierr = VecSetValues( *vecAbove, 1, &diagIndAbove, &diag, ADD_VALUES );CHKERRQ(ierr);
	*nonzeros = NCPA::max( nvals, *nonzeros );
	ierr = MatRestoreRow( halfQ, Ji, &nvals, &indices, &values );CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin( *vecBelow );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( *vecBelow );CHKERRQ(ierr);
	ierr = VecAssemblyBegin( *vecAbove );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( *vecAbove );CHKERRQ(ierr);

	ierr = MatDestroy( &halfQ );CHKERRQ(ierr);
	*nonzeros = NZvec / 10;
	return 1;
	

	/* For order 1,1 and 2,2 */
	// std::vector< PetscScalar > numerator_coefficients, denominator_coefficients;
	// Vec ek_below, ek_above, Be_below, Be_above;
	// KSP ksp;
	// Mat B, C, Ctrans, *last_q_powers;
	
	// PetscInt ncoeffs = 3;
	// numerator_coefficients.push_back( 1.0 );
	// numerator_coefficients.push_back( 1.25 );
	// numerator_coefficients.push_back( 5.0/16.0 );
	// denominator_coefficients.push_back( 1.0 );
	// denominator_coefficients.push_back( 0.75 );
	// denominator_coefficients.push_back( 1.0/16.0 );

	// PetscInt ncoeffs = 2;
	// numerator_coefficients.push_back( 1.0 );
	// numerator_coefficients.push_back( 0.75 );
	// denominator_coefficients.push_back( 1.0 );
	// denominator_coefficients.push_back( 0.25 );

	// last_q_powers = NULL;
	// create_matrix_polynomial( ncoeffs, Q, &last_q_powers );
	// generate_polymatrices( last_q_powers, ncoeffs, NZvec, numerator_coefficients, denominator_coefficients, &B, &C );

	// // create index vectors
	// ierr = VecCreate( PETSC_COMM_SELF, &ek_below );CHKERRQ(ierr);
	// ierr = VecSetType( ek_below, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetSizes( ek_below, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( ek_below );CHKERRQ(ierr);
	// ierr = VecSet( ek_below, 0 );CHKERRQ(ierr);
	// ierr = VecCreate( PETSC_COMM_SELF, &ek_above );CHKERRQ(ierr);
	// ierr = VecSetType( ek_above, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetSizes( ek_above, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( ek_above );CHKERRQ(ierr);
	// ierr = VecSet( ek_above, 0 );CHKERRQ(ierr);

	// if (Ji > 0) {
	// 	ierr = VecSetValue( ek_below, Ji-1, 1, INSERT_VALUES );CHKERRQ(ierr);
	// }
	// ierr = VecSetValue( ek_above, Ji, 1, INSERT_VALUES );CHKERRQ(ierr);
	// ierr = VecAssemblyBegin( ek_above );CHKERRQ(ierr);
	// ierr = VecAssemblyEnd( ek_above );CHKERRQ(ierr);
	// ierr = VecAssemblyBegin( ek_below );CHKERRQ(ierr);
	// ierr = VecAssemblyEnd( ek_below );CHKERRQ(ierr);

	// ierr = VecCreate( PETSC_COMM_SELF, &Be_below );CHKERRQ(ierr);
	// ierr = VecSetType( Be_below, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetSizes( Be_below, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( Be_below );CHKERRQ(ierr);
	// ierr = VecAssemblyBegin( Be_below );CHKERRQ(ierr);
	// ierr = VecAssemblyEnd( Be_below );CHKERRQ(ierr);

	// ierr = VecCreate( PETSC_COMM_SELF, &Be_above );CHKERRQ(ierr);
	// ierr = VecSetType( Be_above, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetSizes( Be_above, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( Be_above );CHKERRQ(ierr);
	// ierr = VecAssemblyBegin( Be_above );CHKERRQ(ierr);
	// ierr = VecAssemblyEnd( Be_above );CHKERRQ(ierr);

	// // setup RHS and result vectors
	// //ierr = MatCreateTranspose( B, &Btrans );CHKERRQ(ierr);
	// ierr = MatCreateTranspose( C, &Ctrans );CHKERRQ(ierr);
	// ierr = MatMultTranspose( B, ek_below, Be_below );CHKERRQ(ierr);
	// ierr = MatMultTranspose( B, ek_above, Be_above );CHKERRQ(ierr);
	// ierr = VecCreate( PETSC_COMM_SELF, vecBelow );CHKERRQ(ierr);
	// ierr = VecSetType( *vecBelow, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( *vecBelow );CHKERRQ(ierr);
	// ierr = VecSetSizes( *vecBelow, PETSC_DECIDE, NZvec );CHKERRQ(ierr);
	// ierr = VecCreate( PETSC_COMM_SELF, vecAbove );CHKERRQ(ierr);
	// ierr = VecSetType( *vecAbove, VECSEQ );CHKERRQ(ierr);
	// ierr = VecSetFromOptions( *vecAbove );CHKERRQ(ierr);
	// ierr = VecSetSizes( *vecAbove, PETSC_DECIDE, NZvec );CHKERRQ(ierr);

	// // Set up solution
	// ierr = KSPCreate( PETSC_COMM_SELF, &ksp );CHKERRQ(ierr);
	// ierr = KSPSetOperators( ksp, Ctrans, Ctrans );CHKERRQ(ierr);
	// ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
	// ierr = KSPSolve( ksp, Be_below, *vecBelow );CHKERRQ(ierr);
	// ierr = KSPSetOperators( ksp, Ctrans, Ctrans );CHKERRQ(ierr);
	// ierr = KSPSolve( ksp, Be_above, *vecAbove );CHKERRQ(ierr);

	// // clean up
	// ierr = KSPDestroy( &ksp );CHKERRQ(ierr);
	// ierr = VecDestroy( &ek_below );CHKERRQ(ierr);
	// ierr = VecDestroy( &ek_above );CHKERRQ(ierr);
	// ierr = VecDestroy( &Be_above );CHKERRQ(ierr);
	// ierr = VecDestroy( &Be_below );CHKERRQ(ierr);
	// // ierr = MatDestroy( &Btrans );CHKERRQ(ierr);
	// ierr = MatDestroy( &Ctrans );CHKERRQ(ierr);
	// ierr = MatDestroy( &B );CHKERRQ(ierr);
	// ierr = MatDestroy( &C );CHKERRQ(ierr);
	// delete_matrix_polynomial( ncoeffs, &last_q_powers );
	// *nonzeros = NZvec;
	// return 1;

}


void NCPA::EPadeSolver::absorption_layer( double lambda, double *z, int NZ, double *layer ) {
	double thickness = NCPA::min<double>(
		NCPAPROP_EPADE_PE_ABSORBING_LAYER_MAX_THICKNESS_METERS,
		lambda * NCPAPROP_EPADE_PE_ABSORBING_LAYER_WAVELENGTH_MULTIPLIER );
	double z_t = z[NZ-1] - thickness;
	for (int i = 0; i < NZ; i++) {
		layer[ i ] = absorption_layer_mu * std::exp( (z[i]-z_t) / thickness );
	}
}

// @todo make static
int NCPA::EPadeSolver::get_starter_gaussian( size_t NZ, double *z,
	double zs, double k0, int ground_index, Vec *psi ) {

	double fac = 2.0;
	//double kfac = k0 / fac;
	PetscScalar tempval;
	PetscErrorCode ierr;

	ierr = VecCreate( PETSC_COMM_SELF, psi );CHKERRQ(ierr);
	ierr = VecSetSizes( *psi, PETSC_DECIDE, NZ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( *psi ); CHKERRQ(ierr);
	ierr = VecSet( *psi, 0.0 );

	for (PetscInt i = 0; i < (PetscInt)NZ; i++) {
		//if (z[i] >= zg) {
			tempval = -( k0*k0/fac/fac ) * (z[i] - zs) * (z[i] - zs);
			tempval = sqrt( k0/fac ) * exp( tempval );
			ierr = VecSetValues( *psi, 1, &i, &tempval, INSERT_VALUES );CHKERRQ(ierr);
		//}
	}
	ierr = VecAssemblyBegin( *psi );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( *psi );CHKERRQ(ierr);
	return 1;
}

// @todo make static
void NCPA::EPadeSolver::make_point_source( size_t NZ, double *z, double zs,
		double z_ground, std::complex<double> *source ) {
//	std::memset( source, 0, NZ * sizeof(std::complex<double>) );
	std::fill( source, source+NZ, std::complex<double>{} );
	size_t nzsrc = NCPA::find_closest_index<double>( z, NZ, zs );
	while (z[nzsrc] < z_ground) {
		nzsrc++;
	}
	source[ nzsrc ].real( 1.0 );
}

void NCPA::EPadeSolver::read_line_source_from_file( size_t NZ, double *z,
	double z_ground, const std::string &filename,
	std::complex<double> *source ) {

	std::vector< std::string >::const_iterator cit;

	// get the lines
	std::string delimiters = ":,= ";
	std::string headerchars = "#";
	std::vector< std::vector< std::string > > contents;
	std::vector< std::string > headerlines;
	NCPA::read_text_columns_from_file_with_header(
		filename, contents, headerlines, delimiters, headerchars );

	// first, parse the header for units information
	NCPA::units_t file_z_units = NCPAPROP_EPADE_PE_UNITS_Z;
	for (cit = headerlines.cbegin(); cit != headerlines.cend(); ++cit) {
		std::string thisline = *cit;
		if (thisline.find( "#%" ) == 0) {
			thisline.erase( 0, 2 );
			thisline = NCPA::deblank( thisline );
			std::vector<std::string> fields = NCPA::split(
				thisline, delimiters );
			if (fields.size() < 2) {
				std::cerr << "Line source file descriptive header line "
						  << thisline << " has no delimiter characters ("
						  << delimiters << "), ignoring" << std::endl;
			} else {
				units_t tempunits = NCPA::Units::fromString( fields[1] );
				if (tempunits == UNITS_NONE) {
					std::cerr << "Unrecognized units " << fields[1]
							  << ", ignoring" << std::endl;
				} else {
					switch ((fields[0])[0]) {
						case 'z':
						case 'Z':
							file_z_units = tempunits;
							break;
						default:
							std::cerr << "Unrecognized parameter tag "
									<< (fields[0])[0]
									<< ", must be in [Zz].  Ignoring"
									<< std::endl;
					}
				}
			}
		}
	}

	// now get the column contents
	size_t ncols = contents.size();
	size_t nvals = contents[ 0 ].size();
	bool complex_in = (ncols == 3);

	double *z_orig = NCPA::zeros<double>( nvals ),
		   *r_orig = NCPA::zeros<double>( nvals ),
		   *i_orig = NCPA::zeros<double>( nvals );
	for (size_t ii = 0; ii < nvals; ii++) {
		z_orig[ ii ] = std::stod( contents[ 0 ][ ii ] );
		r_orig[ ii ] = std::stod( contents[ 1 ][ ii ] );
		if (complex_in) {
			i_orig[ ii ] = std::stod( contents[ 2 ][ ii ] );
		}
	}

//	std::memset( source, 0, NZ*sizeof(std::complex<double>) );
	// convert units
	NCPA::Units::convert(z_orig,nvals,file_z_units,NCPAPROP_EPADE_PE_UNITS_Z,z_orig);
	std::fill( source, source+NZ, std::complex<double>{} );
	interpolate_complex( nvals, z_orig, r_orig, i_orig, NZ, z, source );
	delete [] z_orig;
	delete [] r_orig;
	delete [] i_orig;

}


int NCPA::EPadeSolver::get_starter_self( size_t NZ, double *z,
	std::complex<double> *source, double k0, Mat *qpowers, size_t npade,
	Vec *psi ) {

	Vec rhs, ksi, Bksi, tempvec;
	Mat /*A, AA,*/ B, C;
	KSP /*ksp,*/ ksp2;
	PetscScalar I( 0.0, 1.0 ), tempsc, zeroval = 0.0;
	// PetscInt ii, Istart, Iend;
	PetscErrorCode ierr;
	
	// create rhs vector
	ierr = VecCreate( PETSC_COMM_SELF, &rhs );CHKERRQ(ierr);
	ierr = VecSetSizes( rhs, PETSC_DECIDE, NZ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( rhs );CHKERRQ(ierr);
	ierr = VecSet( rhs, 0.0 );CHKERRQ(ierr);
	
	// find closest index to zs. Make sure the picked point is above
	// the ground surface if we're working in absolute elevation.  If
	// we're in relative elevation, the ground is at 0 by definition
	// size_t nzsrc;
	// if (absolute) {
	// 	nzsrc = NCPA::find_closest_index<double>( z, NZ, zs+z_ground );
	// 	while (z[nzsrc] < z_ground) {
	// 		nzsrc++;
	// 	}
	// } else {
	// 	nzsrc = NCPA::find_closest_index<double>( z, NZ, zs );
	// }

	double h = z[1] - z[0];
	PetscScalar hinv = 1.0 / h;
	PetscInt *indices = NCPA::index_vector<PetscInt>( NZ );
	// PetscInt ps_nzsrc = nzsrc;
	// ierr = VecSetValues( rhs, 1, &ps_nzsrc, &hinv, INSERT_VALUES );CHKERRQ(ierr);
	ierr = VecSetValues( rhs, NZ, indices, source, INSERT_VALUES );CHKERRQ(ierr);
	ierr = VecScale( rhs, hinv );CHKERRQ(ierr);

	// solve first part (Eq. 26)
	ierr = VecDuplicate( rhs, &ksi );CHKERRQ(ierr);
	ierr = VecCopy( rhs, ksi );CHKERRQ(ierr);
	
	// get starter
	std::cout << "Finding ePade starter coefficients..." << std::endl;
	double r_ref = 2 * PI / k0;
	std::vector< PetscScalar > P, Q;
	std::vector< PetscScalar > taylor1 = 
		taylor_sqrt_1pQ_exp_id_sqrt_1pQ_m1( 2*npade, k0*r_ref );
	calculate_pade_coefficients( &taylor1, npade, npade+1, &P, &Q );

	generate_polymatrices( qpowers, npade, NZ, P, Q, &B, &C );
	PetscScalar hank_inv = pow( sqrt( 2.0 / ( PI * k0 * r_ref ) ) * exp( I * (k0 * r_ref - PI/4.0 ) ),
		-1.0 );

	// Original Matlab: psi = AA * ( C \ (B * ksi) ) / hank
	// compute product of B and ksi
	ierr = VecDuplicate( ksi, &Bksi );
	ierr = VecDuplicate( ksi, &tempvec );
	ierr = VecDuplicate( ksi, psi );
	ierr = MatMult( B, ksi, Bksi );

	// solve for tempvec = C \ Bksi
	ierr = KSPCreate( PETSC_COMM_WORLD, &ksp2 );CHKERRQ(ierr);
	ierr = KSPSetOperators( ksp2, C, C );CHKERRQ(ierr);
	ierr = KSPSetFromOptions( ksp2 );CHKERRQ(ierr);
	ierr = KSPSolve( ksp2, Bksi, *psi );CHKERRQ(ierr);
	
	// multiply and scale
	ierr = VecScale( *psi, hank_inv );CHKERRQ(ierr);


	// clean up
	ierr = MatDestroy( &B );CHKERRQ(ierr);
	ierr = MatDestroy( &C );CHKERRQ(ierr);
	ierr = VecDestroy( &rhs );CHKERRQ(ierr);
	ierr = VecDestroy( &ksi );CHKERRQ(ierr);
	ierr = VecDestroy( &Bksi );CHKERRQ(ierr);
	ierr = VecDestroy( &tempvec );CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp2 );CHKERRQ(ierr);

	return 1;
}






int NCPA::EPadeSolver::calculate_pade_coefficients( std::vector<PetscScalar> *c, 
	int n_numerator, int n_denominator, std::vector<PetscScalar> *numerator_coefficients,
	std::vector<PetscScalar> *denominator_coefficients ) {

	// sanity checks
	if (n_denominator < n_numerator) {
		std::cerr << "Denominator count must be >= numerator count for Pade calculation" << std::endl;
		exit(0);
	}
	int n = n_numerator - 1;    // numerator order
	int m = n_denominator - 1;  // denominator order
	int N = n + m;
	int n_taylor = c->size();
	if (n_taylor < (N+1)) {
		std::cerr << "Count of Taylor series must be at least " << (N+1) << " for numerator count "
				  << n_numerator << " and denominator count " << n_denominator << std::endl;
		exit(0);
	}

	//double delta = k0 * dr;
	std::complex<double> j( 0.0, 1.0 );
	PetscErrorCode ierr;
	PetscInt Istart, Iend, ii, jj, *indices;
	// PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar tempsc, *contents;
	Mat A;
	Vec x, y;
	KSP ksp;

	// Create and populate matrix system
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, N, N, n_denominator, NULL, &A );CHKERRQ(ierr);
	ierr = MatSetFromOptions( A );CHKERRQ(ierr);
	ierr = MatZeroEntries( A );CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
	tempsc = -1.0;
    for (ii = Istart; ii < min(n,Iend); ii++) {
    	ierr = MatSetValues(A,1,&ii,1,&ii,&tempsc,INSERT_VALUES);CHKERRQ(ierr);
    }
    for (ii = Istart; ii < Iend; ii++) {
    	for (jj = n; jj <= min(Iend-1,ii+n); jj++) {
    		tempsc = c->at( ii-jj+n );
    		ierr = MatSetValues(A,1,&ii,1,&jj,&tempsc,INSERT_VALUES);CHKERRQ(ierr);
    	} 
    }


    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    // outputSparseMat( A, N, "pade_matrix.dat" );

	// setup right-side vector
	ierr = VecCreate( PETSC_COMM_SELF, &x );CHKERRQ(ierr);
	ierr = VecSetSizes( x, PETSC_DECIDE, N );CHKERRQ(ierr);
	ierr = VecSetFromOptions( x ); CHKERRQ(ierr);
	ierr = VecDuplicate( x, &y );CHKERRQ(ierr);

	indices = NCPA::zeros<PetscInt>( N );
	for (ii = 0; ii < N; ii++) {
		tempsc = -c->at( ii+1 );
		// std::cout << "y[ " << ii << " ] = -c[ " << ii+1 << " ] = " << tempsc << std::endl;
		ierr = VecSetValues( y, 1, &ii, &tempsc, INSERT_VALUES );CHKERRQ(ierr);
		indices[ ii ] = ii;
	}
	ierr = VecAssemblyBegin( y );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( y );CHKERRQ(ierr);

	ierr = VecSet( x, 0.0 );CHKERRQ(ierr);

	// solve
	ierr = KSPCreate( PETSC_COMM_WORLD, &ksp );CHKERRQ(ierr);
	ierr = KSPSetOperators( ksp, A, A );CHKERRQ(ierr);
	ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
	ierr = KSPSolve( ksp, y, x );CHKERRQ(ierr);

	// populate P and Q vectors. Q is denominator coefficients (b), P is numerator coefficients (a)
	contents = NCPA::zeros<PetscScalar>( N );
	ierr = VecGetValues( x, N, indices, contents );

	numerator_coefficients->clear();
	numerator_coefficients->push_back( c->at(0) );
	for (ii = 0; ii < n; ii++) {
		numerator_coefficients->push_back( contents[ ii ] );
	}
	denominator_coefficients->clear();
	denominator_coefficients->push_back( 1.0 );
	for (ii = n; ii < N; ii++) {
		denominator_coefficients->push_back( contents[ ii ] );
	}
	delete [] contents;
	delete [] indices;
	
	// clean up memory
	ierr = KSPDestroy( &ksp );CHKERRQ(ierr);
	ierr = VecDestroy( &x );CHKERRQ(ierr);
	ierr = VecDestroy( &y );CHKERRQ(ierr);
	ierr = MatDestroy( &A );CHKERRQ(ierr);
	return 0;
}

// Uses the recurrence relation derived by Assink to modify the Taylor series for
// F=exp[ i*d*( sqrt(1+Q) - 1 ) ] to F=exp[ i*d*( sqrt(1+Q) - 1 ) ] / sqrt( 1+Q )
std::vector<PetscScalar> NCPA::EPadeSolver::taylor_sqrt_1pQ_exp_id_sqrt_1pQ_m1( 
	int N, double delta ) {

	std::complex<double> j( 0.0, 1.0 );

	// first get the original series
	std::vector<PetscScalar> c = taylor_exp_id_sqrt_1pQ_m1( N, delta );

	// Now modify
	std::vector<PetscScalar> d( N, 1.0 );
	for (int m = 1; m < N; m++) {
		d[ m ] = (j*delta*c[m-1] - ( ((double)(1 + 2*(m-1))) * d[m-1])) / (2.0*m);
	}

	return d;
}

// Uses the recurrence relation in Roberts & Thompson (2013, eq. 16) to compute the
// Taylor series coefficients for F=exp[ i*d*( sqrt(1+Q) - 1 ) ] up to order N-1 
// (i.e. the first N terms)
std::vector<PetscScalar> NCPA::EPadeSolver::taylor_exp_id_sqrt_1pQ_m1( int N, double delta ) {
	std::complex<double> j( 0.0, 1.0 );
	std::vector<PetscScalar> c( N, 1.0 );
	//std::memset( c, 0, N*sizeof(PetscScalar) );
	//c[ 0 ] = 1.0;
	c[ 1 ] = j * 0.5 * delta;
	for ( int idx= 2; idx < N; idx++) {
		double dm = (double)(idx - 1);
		c[ idx ] = -((2.0*dm - 1.0) / (2.0*dm + 2.0)) * c[idx-1]
				   - (delta*delta / (4.0*dm*(dm+1.0))) * c[idx-2];
	}

	return c;
}

// Uses a calculated recursion relation to compute the Taylor series coefficients for
// G=(1+Q)^-0.25 up to order N-1 (i.e. the first N terms)
std::vector<PetscScalar> NCPA::EPadeSolver::taylor_1pQ_n025( int N ) {
	std::vector<PetscScalar> c( N, 1.0 );
	c[ 1 ] = -0.25;
	for (int idx = 2; idx < N; idx++) {
		double dn = double(idx);
		c[ idx ] = -((4*dn)-3) * c[ idx-1 ] / (4*dn); 
	}

	return c;
}

// Uses a calculated recursion relation to compute the Taylor series coefficients for
// G=(1+Q)^0.25 up to order N-1 (i.e. the first N terms)
std::vector<PetscScalar> NCPA::EPadeSolver::taylor_1pQ_025( int N ) {
	std::vector<PetscScalar> c( N, 1.0 );
	c[ 1 ] = 0.25;
	for (int idx = 1; idx < N; idx++) {
		double dn = double(idx);
		c[ idx ] = c[ idx-1 ] * (1.0/dn) * -((4.0*dn - 5.0)/4.0);
	}

	return c;
}

std::vector<PetscScalar> NCPA::EPadeSolver::taylor_1pQpid_n025( int N, double delta ) {
	std::vector<PetscScalar> c( N, 0.0 );
	std::complex< double > J( 0.0, 1.0 );
	c[ 0 ] = std::pow( 1.0 + J*delta, -0.25 );
	for (int idx = 1; idx < N; idx++) {
		double dn = double(idx);
		c[ idx ] = c[ idx-1 ] * (1.0/dn) * (-(4.0*(dn-1.0)+1.0) / 4) * std::pow( 1.0 + J*delta, -1.0 );
	}

	return c;
}



void NCPA::EPadeSolver::output1DTL( std::string filename, double calc_az, bool append ) {
	std::ofstream out_1d;
	if (append) {
		out_1d.open( filename, std::ofstream::out | std::ofstream::app );
		out_1d << std::endl;
	} else {
		out_1d.open( filename, std::ofstream::out | std::ofstream::trunc );
	}
	for (size_t i = 0; i < (NR-1); i++) {
		out_1d << r[ i ]/1000.0 << " "
			   << calc_az << " "
//			   << z[ zgi_r[ i ] ] << " "
			   << tl[ zgi_r[ i ] ][ i ].real() << " "
			   << tl[ zgi_r[ i ] ][ i ].imag() << std::endl;
	}
	out_1d.close();
}

void NCPA::EPadeSolver::output2DTL( std::string filename ) {
	std::ofstream out_2d( filename, std::ofstream::out | std::ofstream::trunc );
	for (size_t i = 0; i < (NR-1); i++) {
		for (size_t j = 0; j < NZ; j += NCPAPROP_EPADE_PE_2D_OUTPUT_Z_STEP) {
			out_2d  << r[ i ]/1000.0 << " "
					<< z[ j ]/1000.0 << " "
					<< tl[ j ][ i ].real() << " "
					<< tl[ j ][ i ].imag() << std::endl;
		}
		out_2d << std::endl;
	}
	out_2d.close();
}

void NCPA::EPadeSolver::set_1d_output( bool tf ) {
	write1d = tf;
}

/*
Broadband internal header format:

uint32_t n_az
uint32_t n_f
uint32_t precision_factor
int64_t  az[ 0 ] * precision_factor
  ...
int64_t  az[ n_az-1 ] * precision_factor
int64_t  f[ 0 ] * precision_factor
  ...
int64_t  f[ n_f-1 ] * precision_factor
[ body ]
*/
void NCPA::EPadeSolver::write_broadband_header( std::string filename, double *az_vec, size_t n_az, 
	double *f_vec, size_t n_f, unsigned int precision_factor ) {

	size_t i = 0;

	// open the file, truncating it if it exists
	std::ofstream ofs( filename, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary );
	if (!ofs.good()) {
		throw std::runtime_error( "Error opening file to initialize:" + filename );
	}

	size_t buf_size = n_az;
	if (n_f > buf_size) {
		buf_size = n_f;
	}
	int64_t *buffer = NCPA::zeros<int64_t>( buf_size );

	// write header starting with vector sizes and multiplicative factor
	uint32_t holder = n_az;
	ofs.write( (char*)(&holder), sizeof( uint32_t ) );
	holder = n_f;
	ofs.write( (char*)(&holder), sizeof( uint32_t ) );
	holder = precision_factor;
	ofs.write( (char*)(&holder), sizeof( uint32_t ) );

	for (i = 0; i < n_az; i++) {
		buffer[ i ] = (int64_t)std::lround( az_vec[ i ] * (double)precision_factor );
	}
	ofs.write( (char*)buffer, n_az * sizeof( int64_t ) );
	std::memset( buffer, 0, buf_size * sizeof( int64_t ) );
	for (i = 0; i < n_f; i++) {
		buffer[ i ] = (int64_t)std::lround( f_vec[ i ] * (double)precision_factor );
	}
	ofs.write( (char*)buffer, n_f * sizeof( int64_t ) );
	ofs.close();

	delete [] buffer;
}

/*
Broadband body format:
foreach (az)
  foreach (freq)
    int64_t  az                       * precision_factor
    int64_t  freq                     * precision_factor
    uint32_t n_z
    uint32_t n_range
    int64_t  z[ 0 ]                   * precision_factor
      ...
    int64_t  z[ n_z-1 ]               * precision_factor
    int64_t  range[ 0 ]               * precision_factor
      ...
    int64_t  range[ n_range-1 ]       * precision_factor
    int64_t  Re{ TL[ z[0] ][ r[0] ] } * precision_factor
	int64_t  Im{ TL[ z[0] ][ r[0] ] } * precision_factor
	int64_t  Re{ TL[ z[0] ][ r[1] ] } * precision_factor
	int64_t  Im{ TL[ z[0] ][ r[1] ] } * precision_factor
	  ...
	int64_t  Re{ TL[ z[0] ][ r[n_range-1] ] } * precision_factor
	int64_t  Im{ TL[ z[0] ][ r[n_range-1] ] } * precision_factor
	int64_t  Re{ TL[ z[1] ][ r[0] ] } * precision_factor
	int64_t  Im{ TL[ z[1] ][ r[0] ] } * precision_factor
	int64_t  Re{ TL[ z[1] ][ r[1] ] } * precision_factor
	int64_t  Im{ TL[ z[1] ][ r[1] ] } * precision_factor
	  ...
*/
void NCPA::EPadeSolver::write_broadband_results( std::string filename, double this_az, double this_f, 
	double *r_vec, size_t n_r, double *z_vec, size_t n_z, std::complex< double > **tloss_mat, 
	unsigned int precision_factor ) {

	n_r--;    // last range step is invalid

	std::ofstream ofs( filename, std::ofstream::out | std::ofstream::app | std::ofstream::binary );
	if (!ofs.good()) {
		throw std::runtime_error( "Error opening file to append: " + filename );
	}

	// write az, freq, n_z, n_range
	int64_t holder = (int64_t)std::lround( this_az * (double)precision_factor );
	ofs.write( (char*)(&holder), sizeof( int64_t ) );
	holder = (int64_t)std::lround( this_f * (double)precision_factor );
	ofs.write( (char*)(&holder), sizeof( int64_t ) );
	uint32_t uholder = (uint32_t)n_z;
	ofs.write( (char*)(&uholder), sizeof( uint32_t ) );
	uholder = (uint32_t)n_r;
	ofs.write( (char*)(&uholder), sizeof( uint32_t ) );

	// z and r sizes and vectors
	size_t buf_size = n_r;
	if (n_z > buf_size) {
		buf_size = n_z;
	}
	int64_t *buffer = NCPA::zeros<int64_t>( buf_size );
	size_t i, j;
	for (i = 0; i < n_z; i++) {
		buffer[ i ] = (int64_t)std::lround( z_vec[ i ] * (double)precision_factor );
	}
	ofs.write( (char*)buffer, n_z * sizeof( int64_t ) );
	for (i = 0; i < n_r; i++) {
		buffer[ i ] = (int64_t)std::lround( r_vec[ i ] * (double)precision_factor );
	}
	ofs.write( (char*)buffer, n_r * sizeof( int64_t ) );
	for (i = 0; i < n_z; i++) {
		for (j = 0; j < n_r; j++) {
			holder = (int64_t)std::lround( tloss_mat[ i ][ j ].real() * (double)precision_factor );
			ofs.write( (char *)(&holder), sizeof( int64_t ) );
			holder = (int64_t)std::lround( tloss_mat[ i ][ j ].imag() * (double)precision_factor );
			ofs.write( (char *)(&holder), sizeof( int64_t ) );
		}
	}
	ofs.close();
	delete [] buffer;
}

int NCPA::EPadeSolver::zero_below_ground( Mat *q, int NZ, PetscInt ground_index ) {

	PetscErrorCode ierr;
	PetscInt nNonZero, ii, jj, kk;
	const PetscInt *indices;
	//const PetscScalar *values;
	PetscScalar zero = 0.0;

	std::vector< PetscInt > rows, cols;

	for (ii = 0; ii < ground_index; ii++) {
		ierr = MatGetRow( *q, ii, &nNonZero, &indices, PETSC_NULL );CHKERRQ(ierr);
		for (jj = 0; jj < nNonZero; jj++) {
			rows.push_back( ii );
			cols.push_back( indices[ jj ] );
		}
		ierr = MatRestoreRow( *q, ii, &nNonZero, &indices, PETSC_NULL );CHKERRQ(ierr);
	}
	if (ground_index > 0) {
		rows.push_back( ground_index );
		cols.push_back( ground_index-1 );
	}

	for (ii = 0; ii < (int)(rows.size()); ii++) {
		jj = rows[ ii ];
		kk = cols[ ii ];
		ierr = MatSetValues( *q, 1, &jj, 1, &kk, &zero, INSERT_VALUES );CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin( *q, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);
	ierr = MatAssemblyEnd( *q, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);

	return 1;
}

// @todo make static
void NCPA::EPadeSolver::write_topography( std::string filename,
	double azimuth, double max_range, double dr ) {

	double r;
	std::ofstream outfile( filename, std::ofstream::out | std::ofstream::app );
	if (!outfile.good()) {
		throw std::runtime_error( "Error opening file " + filename + " to write topography" );
	}

	for (r = 0.0; r <= max_range; r += dr) {
		outfile << azimuth << " " << r / 1000.0 << " "
				<< atm_profile_2d->get_interpolated_ground_elevation( r ) 
				<< std::endl;
	}

	outfile.close();
}

void NCPA::EPadeSolver::write_source( const std::string &filename,
		const std::complex<double> *source, const double *z, size_t NZ ) const {
	std::ofstream out(filename,std::ios_base::out);
	for (size_t i = 0; i < NZ; i++) {
		out << z[i] << " " << source[i].real() << " " << source[i].imag() << std::endl;
	}
	out.close();
}

void NCPA::EPadeSolver::calculate_effective_sound_speed(
	NCPA::Atmosphere2D *atm, double az, const std::string &new_key ) {

	// first: was it given explicitly using column "CEFF"?
	if (atm->contains_vector( 0.0, "CEFF" )) {
		atm->convert_property_units( "CEFF", NCPAPROP_EPADE_PE_UNITS_C );
		atm->copy_vector_property( "CEFF", new_key );

	// do we have the wind speed and direction?
	} else if (atm->contains_vector( 0.0, "_WS_")
		&& atm->contains_vector( 0.0, "_WD_") ) {
		atm->calculate_wind_component( "_WC_", "_WS_", "_WD_", az );
		atm->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );
		atm->remove_property( "_WC_" );
	} else {
		std::ostringstream oss;
		oss << "Cannot calculate effective sound speed, necessary components not found."
			<< std::endl
			<< "Input atmosphere must have one of:" << std::endl
			<< "  CEFF column, or" << std::endl
			<< "  WS and WD columns for wind speed and direction, or" << std::endl
			<< "  U and V columns for zonal and meridional wind vectors." << std::endl;

		throw std::runtime_error( oss.str() );
	}
}

void NCPA::EPadeSolver::calculate_turbulence( double r,
		size_t nz, double *z, double k_a, size_t ground_index,
		double *&mu ) const {

	size_t i, j, nt;
	nt = turbulence->size();

	gsl_vector_set_zero( t_vec1 );
	gsl_vector_set_zero( t_vec_mu );
	gsl_matrix_set_zero( t_mat1 );

	for (i = 0; i < nt; i++) {
		// vec1->set( 0, i, turbulence->get_G( i ) );
		gsl_vector_set( t_vec1, i, turbulence->get_G( i ) );
		for (j = 0; j < nz; j++) {
			double temp = r * turbulence->get_k( i ).real()
						+ turbulence->get_alpha( i )
						+ turbulence->get_k( i ).imag() * z[ j ];
			gsl_matrix_set( t_mat1, j, i, std::cos( temp ) );
			// mat1->set( i, j, std::cos( temp ) );
			// ofs << i << " " << j << " " << mat1->get( i, j ) << std::endl;
		}
	}

	gsl_blas_dgemv( CblasNoTrans, 1.0, t_mat1, t_vec1, 0.0, t_vec_mu );
	for (j = 0; j < nz; j++) {
		if (j >= ground_index) {
			mu[ j ] = gsl_vector_get( t_vec_mu, j );
		} else {
			mu[ j ] = 0.0;
		}
	}
}

// int NCPA::EPadeSolver::calculate_turbulence( double r,
// 		size_t nz, double *z, double k_a,
// 		double *&mu ) const {

// 	// Calculates Eq. J.24 from Salomons.
// 	PetscInt i, j, nt;
// 	nt = turbulence->size();

// 	Mat mat1;
// 	Vec vec_mu, vec1;
// 	PetscErrorCode ierr;

// 	// set up vector and matrix objects
// 	ierr = MatCreateSeqDense( PETSC_COMM_SELF, nz, nt, NULL, &mat1 );CHKERRQ(ierr);
// 	ierr = MatSetFromOptions( mat1 );CHKERRQ(ierr);
// 	ierr = VecCreate( PETSC_COMM_SELF, &vec1 );CHKERRQ(ierr);
// 	ierr = VecSetSizes( vec1, PETSC_DECIDE, nt );CHKERRQ(ierr);
// 	ierr = VecSetFromOptions( vec1 );CHKERRQ(ierr);
// 	ierr = VecSet( vec1, 0.0 );CHKERRQ(ierr);


// 	ierr = VecCreate( PETSC_COMM_SELF, &vec_mu );CHKERRQ(ierr);
// 	ierr = VecSetSizes( vec_mu, PETSC_DECIDE, nz );CHKERRQ(ierr);
// 	ierr = VecSetFromOptions( vec_mu );CHKERRQ(ierr);
// 	ierr = VecSet( vec_mu, 0.0 );CHKERRQ(ierr);

// 	// fill vector and matrix
// 	for (i = 0; i < nt; i++) {
// 		PetscScalar tmpG = turbulence->get_G( i );
// 		ierr = VecSetValues( vec1, 1, &i, &tmpG, INSERT_VALUES );CHKERRQ(ierr);
// 		for (j = 0; j < (PetscInt)nz; j++) {
// 			std::complex<double> temp( std::cos(
// 						r * turbulence->get_k( i ).real()
// 						+ turbulence->get_alpha( i )
// 						+ turbulence->get_k( i ).imag() * z[ j ] ),
// 						0.0 );
// 			ierr = MatSetValues( mat1, 1, &j, 1, &i, &temp, INSERT_VALUES );CHKERRQ(ierr);
// 		}
// 	}
// 	// outputVec( vec1, NULL, nt, "vec1_new.dat" );
// 	ierr = MatAssemblyBegin( mat1, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);
// 	ierr = MatAssemblyEnd( mat1, MAT_FINAL_ASSEMBLY );CHKERRQ(ierr);
// 	// if (first_time) {
// 		ierr = VecAssemblyBegin( vec1 );CHKERRQ(ierr);
// 		ierr = VecAssemblyEnd( vec1 );CHKERRQ(ierr);
// 	// }
// 	ierr = VecAssemblyBegin( vec_mu );CHKERRQ(ierr);
// 	ierr = VecAssemblyEnd( vec_mu );CHKERRQ(ierr);

// 	// mat_mu = vec1->multiply( mat1 );
// 	ierr = MatMult( mat1, vec1, vec_mu );CHKERRQ(ierr);
// 	PetscInt *indices = NCPA::index_vector<PetscInt>( nz );
// 	std::complex<double> *cmu = NCPA::zeros<std::complex<double>>( nz );
// 	ierr = VecGetValues( vec_mu, nz, indices, cmu );CHKERRQ(ierr);

// 	for (j = 0; j < (PetscInt)nz; j++) {
// 		mu[ j ] = cmu[ j ].real();
// 	}

// 	delete [] cmu;
// 	ierr = VecDestroy( &vec_mu );CHKERRQ(ierr);
// 	ierr = VecDestroy( &vec1 );CHKERRQ(ierr);
// 	ierr = MatDestroy( &mat1 );CHKERRQ(ierr);

// 	return 1;
// }

void NCPA::EPadeSolver::setup_turbulence(std::vector<double> &rand1,
		std::vector<double> &rand2 ) {
	if (random_turbulence) {
		rand1 = NCPA::random_numbers( turbulence_size );
		rand2 = NCPA::random_numbers( turbulence_size );
	} // otherwise they're already precalculated
	turbulence = new NCPA::Turbulence( turbulence_size );
	turbulence->set_turbulence_scale( Lt );
	// turbulence->set_reference_temperature( T0 );
	turbulence->set_temperature_factor( temperature_factor );
	turbulence->set_velocity_factor( velocity_factor );
	turbulence->set_wavenumbers_log( turbulence_k1,
		turbulence_k2 );
	turbulence->compute_phases( rand1 );
	turbulence->compute();
	turbulence->set_alpha( rand2 );

	t_vec1   = gsl_vector_alloc( turbulence_size );
	t_vec_mu = gsl_vector_alloc( NZ );
	t_mat1   = gsl_matrix_alloc( NZ, turbulence_size );



	// if (turbulence_vec1 != PETSC_NULL) {
	// 	PetscErrorCode ierr = VecDestroy( turbulence_vec1 );CHKERRQ(ierr);
	// 	turbulence_vec1 = PETSC_NULL;
	// }

}

void NCPA::EPadeSolver::cleanup_turbulence() {
	delete turbulence;

	gsl_vector_free( t_vec1 );
	gsl_vector_free( t_vec_mu );
	gsl_matrix_free( t_mat1 );
}
