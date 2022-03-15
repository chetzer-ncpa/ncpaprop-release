#include "GFPESolver.h"
#include <fstream>
#include <cassert>
#include <complex>
#include <stdexcept>
#include <vector>
// #include <map>
#include <set>
#include <utility>
#include "util.h"
#include "ncpa_fft.h"
#include "ncpaprop_common.h"
#include "ncpaprop_atmosphere.h"
#include "matrix.h"
#include "LANLInterpolation.h"





/**********************************************************************
 * Class TransmissionLossField
 *********************************************************************/
NCPA::TransmissionLossField::TransmissionLossField(
	double az, double freq,
	size_t nr, double *r_vec, size_t nz, double *z_vec,
	std::complex<double> **TL ) {

	azimuth_ = az;
	frequency_ = freq;
	r_vec_.reserve( nr );
	z_vec_.reserve( nz );
	TL_mat_ = new NCPA::DenseMatrix<std::complex<double>>( nr, nz );

	for (size_t i = 0; i < nr; i++) {
		r_vec_.push_back( r_vec[ i ] );
		for (size_t j = 0; j < nz; j++) {
			if (i == 0) {
				z_vec_.push_back( z_vec[ j ] );
			}

			TL_mat_->set( i, j, TL[ i ][ j ] );
		}
	}
}


NCPA::TransmissionLossField::TransmissionLossField(
			double az, double freq,
			const std::vector<double> &r_vec,
			const std::vector<double> &z_vec,
			const NCPA::DenseMatrix<std::complex<double>> *TL_mat
			) {
	azimuth_ = az;
	frequency_ = freq;
	r_vec_ = r_vec;
	z_vec_ = z_vec;
	TL_mat_ = new NCPA::DenseMatrix<std::complex<double>>( *TL_mat );
}


NCPA::TransmissionLossField::~TransmissionLossField() {
	delete TL_mat_;
	r_vec_.clear();
	z_vec_.clear();
}


// void NCPA::TransmissionLossField::as_arrays(
// 		double &az, double &freq, size_t &nr, double *&r_vec,
// 		size_t &nz, double *&z_vec, std::complex<double> **&TL) {
// 	az = azimuth_;
// 	freq = frequency_;
// 	nr = r_vec_.size();
// 	r_vec = NCPA::zeros<double>( nr );
// 	nz = z_vec_.size();
// 	z_vec = NCPA::zeros<double>( nz );

// 	TL = NCPA::allocate_matrix<std::complex<double>>( nr, nz );
// 	for (size_t i = 0; i < nr; i++) {
// 		r_vec[ i ] = r_vec_[ i ];
// 		for (size_t j = 0; j < nr; j++) {
// 			if (i == 0) {
// 				z_vec[ j ] = z_vec_[ j ];
// 			}
// 			TL[ i ][ j ] = TL_mat_->get( i, j );
// 		}
// 	}
// }


void NCPA::TransmissionLossField::as_vectors(
		double &az, double &freq,
		std::vector<double> &r_vec, std::vector<double> &z_vec,
		NCPA::DenseMatrix<std::complex<double>> *&TL ) const {
	az = azimuth_;
	freq = frequency_;
	r_vec.assign( r_vec_.begin(), r_vec_.end() );
	z_vec.assign( z_vec_.begin(), z_vec_.end() );
	TL = new NCPA::DenseMatrix<std::complex<double>>( *TL_mat_ );
}




//////////////////////////////////////////////////////////////////







NCPA::GFPESolver::GFPESolver() {
	set_default_values();
}

NCPA::GFPESolver::GFPESolver( NCPA::ParameterSet *param ) {
	set_default_values();

	verbose = !(param->wasFound("quiet"));
	nodelay = param->wasFound("no_delay");

	// read from parameter set object
	use_atm_1d = param->wasFound( "atmosfile" );
	use_atm_2d = param->wasFound( "atmosfile2d" );

	if (use_atm_1d) {
		atm = new NCPA::StratifiedAtmosphere2D(
			param->getString( "atmosfile" ),
			param->getString("atmosheaderfile") );
	} else if (use_atm_2d) {
		atm = new NCPA::ProfileSeriesAtmosphere2D(
			param->getString( "atmosfile2d" ),
			param->getString( "atmosheaderfile" ) );
	}

	// atmospheric units
	units_t u_km = Units::fromString("km"), u_m = Units::fromString("m"),
			u_ms = Units::fromString("m/s"), u_k = Units::fromString("K"),
			u_mbar = Units::fromString("mbar");
	atm->convert_range_units( u_m );
	atm->convert_altitude_units( u_m );
	if (atm->contains_vector(0,"T")) {
		atm->convert_property_units( "T", u_k );
	}
	if (atm->contains_vector(0,"P")) {
		atm->convert_property_units( "P", u_mbar );
	}
	if (atm->contains_vector(0,"U")) {
		atm->convert_property_units( "U", u_ms );
	}
	if (atm->contains_vector(0,"V")) {
		atm->convert_property_units( "V", u_ms );
	}
	if (atm->contains_scalar(0,"Z0")) {
		atm->convert_property_units( "Z0", u_m );
		z_ground = new NCPA::ScalarWithUnits( atm->get(0.0,"Z0"), u_m );
	} else {
		z_ground = new NCPA::ScalarWithUnits(
			atm->get_minimum_altitude(0.0), u_m );
		atm->add_property( "Z0", z_ground->get(), u_m );
	}

	// calculate derived quantities
	// double c0 = 0.0;
	for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm->first_profile();
		 it != atm->last_profile(); ++it) {
		if ( (*it)->contains_vector( "C0" ) ) {
			(*it)->convert_property_units( "C0", Units::fromString( "m/s" ) );
			(*it)->copy_vector_property( "C0", "_C0_" );
			// c0 = atm->get( 0.0, "_C0_", z_ground->get() );
		} else {
			if ( (*it)->contains_vector("P") && (*it)->contains_vector("RHO") ) {
				(*it)->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO",
					Units::fromString( "m/s" ) );
				// c0 = atm->get( 0.0, "_C0_", z_ground->get() );
			} else if ( (*it)->contains_vector("T") ) {
				(*it)->calculate_sound_speed_from_temperature( "_C0_", "T",
					Units::fromString( "m/s" ) );
				// c0 = atm->get( 0.0, "_C0_", z_ground->get() );
			} else if ( (*it)->contains_vector( "CEFF" ) ) {
				// c0 = atm->get( 0.0, "CEFF", z_ground->get() );
			} else {
				throw std::runtime_error( "Cannot calculate static sound speed: None of CEFF, C0, T, or (P and RHO) are specified in atmosphere.");
			}
		}
	}

	// wind speed
	if (atm->contains_vector(0,"WS")) {
		atm->convert_property_units("WS", Units::fromString("m/s") );
		atm->copy_vector_property( "WS", "_WS_" );
	} else if (atm->contains_vector(0,"U")
		&& atm->contains_vector(0,"V")) {
		atm->calculate_wind_speed( "_WS_", "U", "V" );
	}

	// wind direction
	if (atm->contains_vector(0,"WD")) {
		atm->convert_property_units("WD",
			NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
		atm->copy_vector_property( "WD", "_WD_" );
	} else if (atm->contains_vector(0,"U")
		&& atm->contains_vector(0,"V")) {
		atm->calculate_wind_direction( "_WD_", "U", "V" );
	}

	// attenuation
	// if ( attnfile.size() > 0 ) {
	// 	atm->read_attenuation_from_file( "_ALPHA_", param->getString( "attnfile" ) );
	// } else {
		// if (  !(atm->contains_vector(0.0, "T")
		// 		&& atm->contains_vector(0.0, "P")
		// 		&& atm->contains_vector(0.0, "RHO") ) ) {
		// 	std::cout << "At least one of T, P, or RHO is absent, switching to lossless propagation"
		// 			  << std::endl;
		// 	// lossless = true;
		// }
	// }

	// relative humidity
	double Hr = param->getFloat( "humidity" );
	if ( Hr >= 0.0 ) {
		double *Hvec = NCPA::single_valued_vector<double>( atm->nz( 0.0 ), Hr );
		atm->add_property( "_H_", atm->nz( 0.0 ), Hvec, NCPA::UNITS_NONE );
		delete [] Hvec;
	} else if (atm->contains_scalar( 0.0, "H" ) ) {
		Hr = atm->get( 0.0, "H" );
		double *Hvec = NCPA::single_valued_vector<double>( atm->nz( 0.0 ), Hr );
		atm->add_property( "_H_", atm->nz( 0.0 ), Hvec, NCPA::UNITS_NONE );
		delete [] Hvec;
	} else if (atm->contains_vector( 0.0, "H" )) {
		atm->copy_vector_property( "H", "_H_" );
	}

	// frequency
	f_vec.push_back( param->getFloat( "freq" ) );

	// azimuth
	if (param->wasFound("multiprop")) {

		multiprop = true;
		if (param->wasFound( "write_2d_tloss")) {
			if (verbose) {
				std::cout << "2-D output not defined for multiprop, disabling."
					<< std::endl;
			}
			param->removeParameter( "write_2d_tloss" );
		}
		double az_min = param->getFloat( "azimuth_start" );
		double az_max = param->getFloat( "azimuth_end" );
		double d_az   = param->getFloat( "azimuth_step" );
		if (az_max < az_min) {
			throw std::runtime_error( "azimuth_start cannot be greater than azimuth_end" );
		}
		for (double az = az_min; az <= az_max; az += d_az) {
			az_vec.push_back( az );
		}
	} else {
		az_vec.push_back( param->getFloat( "azimuth" ) );
	}

	// suite mode?
	N_suite = param->getInteger( "n_suite" );
	r_suite.clear();
	z_suite.clear();
	if (N_suite > 1) {
		use_suite_mode = true;
		std::vector<std::string>::const_iterator cit;
		if (param->wasFound("r_suite_km")) {
			std::vector<std::string> rstrs
				= NCPA::split( param->getString( "r_suite_km" ), "," );
			for (cit = rstrs.cbegin(); cit != rstrs.cend(); ++cit) {
				r_suite.push_back(
					NCPA::Units::convert( std::stof( *cit ),
						NCPA::UNITS_DISTANCE_KILOMETERS,
						NCPA::UNITS_DISTANCE_METERS ) );
				std::sort( r_suite.begin(), r_suite.end() );
			}
		}
		if (param->wasFound("z_suite_km")) {
			std::vector<std::string> zstrs
				= NCPA::split( param->getString( "z_suite_km" ), "," );
			for (cit = zstrs.cbegin(); cit != zstrs.cend(); ++cit) {
				z_suite.push_back(
					NCPA::Units::convert( std::stof( *cit ),
						NCPA::UNITS_DISTANCE_KILOMETERS,
						NCPA::UNITS_DISTANCE_METERS ) );
				std::sort( z_suite.begin(), z_suite.end() );
			}
		}
	}

	// TL matrices
	TL_.reserve( az_vec.size() * f_vec.size() );

	// distances and heights
	range = new NCPA::ScalarWithUnits(
		param->getFloat( "maxrange_km" ), u_km );
	range->convert_units( u_m );
	if (use_suite_mode && r_suite.size() == 0) {
		r_suite.push_back( range->get() );
	}
	z_source = new NCPA::ScalarWithUnits(
		param->getFloat( "sourceheight_km" ), u_km );
	z_source->convert_units( u_m );
	if (z_source->get() < z_ground->get()) {
		z_source->set_value( z_ground->get() );
		if (verbose) {
			std::cout << "Adjusting source height to ground height of "
					  << z_source->get() << " "
					  << NCPA::Units::toStr( z_source->get_units() )
					  << std::endl;
		}
	}
	z_receiver = new NCPA::ScalarWithUnits(
		param->getFloat( "receiverheight_km" ), u_km );
	z_receiver->convert_units( u_m );
	if (z_receiver->get() < z_ground->get()) {
		z_receiver->set_value( z_ground->get() );
		if (verbose) {
			std::cout << "Adjusting source height to ground height of "
					  << z_receiver->get() << " "
					  << NCPA::Units::toStr( z_receiver->get_units() )
					  << std::endl;
		}
	}
	if (use_suite_mode && z_suite.size() == 0) {
		z_suite.push_back( z_receiver->get() );
	}

	// create turbulence
	use_turbulence = !(param->wasFound( "no_turbulence" ) );
	turbulence_size = (size_t)(param->getInteger( "n_turbulence" ));
	random_turbulence = !(param->wasFound("turbulence_file"));
	if (!random_turbulence)
		turbulence_file = param->getString("turbulence_file");
	// T0 = param->getFloat( "turbulence_ref_temp" );
	Lt = param->getFloat( "turbulence_scale_m" );
	temperature_factor = param->getFloat( "turbulence_t_factor" );
	velocity_factor    = param->getFloat( "turbulence_v_factor" );

	// Ground impedence
	sigma = param->getFloat( "ground_impedence" );

	// resolution factors
	nz_requested = param->getInteger( "nz" );
	dr_requested = param->getFloat( "dr_m" );
	boundary_thickness_requested = param->getFloat( "boundary_layer_m" );
	surface_thickness_requested = param->getFloat( "surface_layer_m" );

	// wavenumber filter
	f1 = param->getFloat( "k_min" );
	f2 = param->getFloat( "k_max" );

	filetag = param->getString( "filetag" );
}


NCPA::GFPESolver::~GFPESolver() {
	delete atm;
	for (std::vector<NCPA::TransmissionLossField *>::iterator it = TL_.begin();
		it != TL_.end(); ++it) {
		delete *it;
	}
	TL_.clear();
	delete range;
	delete z_source;
	delete z_receiver;
	delete z_ground;
	// delete turbulence;
}


void NCPA::GFPESolver::set_default_values() {
	nz_requested = -1;
	dr_requested = -1.0;
	boundary_thickness_requested = -1.0;
	surface_thickness_requested = -1.0;
	atm = NULL;

	// calculation region defaults
	ablthicknessfactor = 50.0;			// ABL thickness factor
	zatmfactor = 100.0;
	deltazfactor = 10.0;
	deltarfactor = 0.1;

	// turbulence defaults
	turbulence_k1 = 0.1;
	turbulence_k2 = 20.0;
	turbulence_size = 20;
	// T0 = 293.0;
	Lt = 100.0;
	temperature_factor = 1.0e-10;
	velocity_factor    = 1.0e-8;
	N_suite = 1;
	filetag = "";

	// flags
	use_turbulence = true;
	random_turbulence = true;
	use_atm_1d = false;
	use_atm_2d = false;
	use_suite_mode = false;
	nodelay = false;
	multiprop = false;
	verbose = true;
}

int NCPA::GFPESolver::solve() {

	double Am = (1.0 - 0.2) / (1000.0 - 30.0);
	double Ab = 1.0 - Am * 1000.0;
	std::vector<double>::const_iterator f, az, vdit;
	std::vector<double>::iterator findit;

	// for suite mode
	std::vector<size_t> r_indices, z_indices;

	std::complex<double> I( 0.0, 1.0 ), one( 1.0, 0.0 );
	size_t i;

	// set up for non-random turbulence if needed
	std::vector<double> rand1, rand2;
	if (use_turbulence && (!random_turbulence)) {
		if (verbose) {
			std::cout << "Reading " << 2*turbulence_size
					  << " values from " << turbulence_file
					  << std::endl;
		}
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

	for (f = f_vec.cbegin(); f != f_vec.cend(); ++f) {

		// ground impedence as a function of frequency
		std::complex<double> Zg = ground_impedence( *f );

		for (az = az_vec.cbegin(); az != az_vec.cend(); ++az) {

			if (verbose) {
				std::cout << "Calculating for frequency " << *f
						  << " Hz, azimuth = " << *az << std::endl;
			}

			// calculate r vector
			std::vector<double> r_vec;

			// average effective sound speed
			calculate_effective_sound_speed( *az, "_CEFF_" );
			double czero = average_value( "_CEFF_", 0.0,
				NCPA::min( z_source->get(), z_receiver->get() ),
				NCPA::max( z_source->get(), z_receiver->get() ) );
			double minzatm = range->get() / std::sqrt( PI * czero );

			// attenuation.  May not have density, so supply humidity as well
			if (atm->contains_vector(0.0, "ALPHA")) {
				atm->copy_vector_property( "ALPHA", "_ALPHA_" );
			} else {
				atm->calculate_attenuation( "_ALPHA_", "T", "P", "RHO",
					"_H_", *f );
			}

			// how much atmosphere do we need?
			double lambda = czero / *f;
			double ablthickness;
			if (boundary_thickness_requested < 0.0) {
				ablthickness = ablthicknessfactor * lambda;
			} else {
				ablthickness = boundary_thickness_requested;
			}

			double zatm;
			if (surface_thickness_requested < 0.0) {
				zatm = NCPA::max<double>( zatmfactor*lambda, minzatm );
			} else {
				zatm = NCPA::max<double>( surface_thickness_requested, minzatm );
			}
			double ztop = zatm + ablthickness;

			// compute new Z vector
			// @todo make sure 2.0*ztop-deltaz is within supplied atmosphere
			size_t Nz;
			double deltaz;
			if (nz_requested < 0) {
				deltaz = lambda / deltazfactor;
				Nz = (size_t)std::exp2(NCPA::nextpow2(ztop/deltaz + 1));
			} else {
				Nz = (size_t)std::exp2(NCPA::nextpow2(nz_requested));
				deltaz = ztop / ((double)Nz);  // approximate
			}
			if (verbose) {
				std::cout << "Using " << Nz << " vertical points."
						  << std::endl;
			}
			size_t Ntrans = 2 * Nz;
			double *z = NCPA::linspace<double>( 0, 2.0*ztop-deltaz, Ntrans );
			deltaz = z[1] - z[0];
			ztop = z[ Nz-1 ];
			double omega = 2.0 * PI * (*f);
			double kzero = omega / czero;
			std::complex<double> betag = kzero / Zg;
			double deltak = (2.0 * PI) / (deltaz * (double)Ntrans);
			double *kprime = NCPA::zeros<double>( Ntrans );

			// find requested Z indices
			z_indices.clear();
			if (use_suite_mode) {
				for (vdit = z_suite.cbegin(); vdit != z_suite.cend(); ++vdit) {
					if (*vdit > ztop) {
						std::cout << "Warning: requested height of "
								  << *vdit << " meters exceeds maximum "
								  << "calculation height of " << ztop
								  << " meters.  Setting requested height to "
								  << ztop << " meters." << std::endl;
					}
					size_t ind = NCPA::find_closest_index( z, Nz, *vdit );
					z_indices.push_back( ind );
				}
			}

			std::complex<double> *vark
				= NCPA::zeros<std::complex<double>>( Ntrans );
			for (i = 0; i < Ntrans; i++) {
				if (i <= Nz) {
					kprime[ i ] = deltak * (double)i;
				} else {
					kprime[ i ] = -deltak * (double)(2*Nz - i);
				}

				// copy this part into r loop
				double alpha =
					boundary_layer( z[ i ], Ab + (*f)*Am,
						ztop, zatm )
					+ atm->get(0.0,"_ALPHA_",z[ i ]);
				double cz = atm->get(0.0,"_CEFF_",z[ i ]);
				std::complex<double> kz( omega / cz, alpha );
				vark[ i ] = kz - kzero;
			}

			double *F = NCPA::zeros<double>( Ntrans );
			wavenumber_filter( Ntrans, kprime, kzero, f1, f2, F );

			// fill r vector
			double dr;
			if (dr_requested < 0.0) {
			 	dr = lambda / deltarfactor;
			} else {
				dr = dr_requested;
			}

			double scalemod = 1.0;
			while (std::log10( dr * std::pow( 10.0, scalemod ) ) < 0.0) {
				scalemod += 1.0;
			}

			// we will round dr to an appropriate number of decimal places
			// by scaling it up, constructing a temporary array of integers,
			// and then scaling them back down
			// @todo use decimal places of requested dr so as not to
			//       lose precision
			double scale = std::pow( 10.0,
				std::floor( std::log10( dr ) ) + scalemod + 1.0 );
			int dr_scaled = (int)std::round( dr * scale );
			int this_r = dr_scaled;
			int rmax_scaled = (int)std::round(range->get() * scale);
			std::vector<int> rscaled;
			while (this_r < rmax_scaled) {
				rscaled.push_back( this_r );
				this_r += dr_scaled;
			}
			rscaled.push_back( rmax_scaled );
			std::vector<int>::iterator iit;

			// If running a suite, insert the exact requested ranges
			if (use_suite_mode) {
				r_indices.clear();
				for (vdit = r_suite.cbegin(); vdit != r_suite.cend(); ++vdit ) {
					rscaled.push_back( (int)std::round( *vdit * scale ) );
				}

				std::sort( rscaled.begin(), rscaled.end() );
				iit = std::unique( rscaled.begin(), rscaled.end() );
				rscaled.resize( std::distance(rscaled.begin(),iit) );

				// find the inserted values after they've been sorted in
				for (vdit = r_suite.cbegin(); vdit != r_suite.cend(); ++vdit ) {
					int val = (int)std::round( *vdit * scale );
					iit = std::find( rscaled.begin(), rscaled.end(), val );
					assert( iit != rscaled.end() );
					r_indices.push_back( std::distance(
						rscaled.begin(), iit ) );
				}
			}

			// de-scale
			r_vec.clear();
			for (iit = rscaled.begin(); iit != rscaled.end(); ++iit) {
				r_vec.push_back( ((double)*iit) / scale );
			}

			for (size_t N = 0; N < N_suite; N++) {
				// set up for loop over range
				std::complex<double> E1, E1F, E2;
				std::complex<double> *E3 = NCPA::zeros<std::complex<double>>( Ntrans );
				std::complex<double> *phi = NCPA::zeros<std::complex<double>>( Ntrans );
				std::complex<double> *phitilde = NCPA::zeros<std::complex<double>>( Ntrans );
				std::complex<double> *phitildeminus = NCPA::zeros<std::complex<double>>( Ntrans );
				std::complex<double> *integrand = NCPA::zeros<std::complex<double>>( Ntrans );
				double *THETA = NCPA::zeros<double>(Ntrans);
				double *dTHETA = NCPA::zeros<double>(Ntrans);
				double *GAMMA = NCPA::zeros<double>(Ntrans);

				// calculate turbulence
				if (use_turbulence) {
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
				}

				NCPA::DenseMatrix< std::complex<double> > *TL_mat
					= new NCPA::DenseMatrix< std::complex<double> >(
						r_vec.size(), Nz );

				std::vector<double>::const_iterator r;
				size_t step_num = 0;
				std::complex<double> scaling = 0.0;
				for (r = r_vec.cbegin(); r != r_vec.cend(); ++r) {

					// first step?
					// std::complex<double> scaling(
					// 	1.0 / std::sqrt( *r ), 0.0 );
					scaling.real( 0.0 );
					scaling.imag( 0.0 );
					if (nodelay) {
						scaling.real( std::pow( *r, -0.5 ) );
					} else {
						std::complex<double> tmpexp1( 0.0, (*r) * kzero );
						std::complex<double> tmpexp2( std::pow( *r, -0.5 ) );
						scaling = std::exp( tmpexp1 ) * tmpexp2;
					}

					if (r == r_vec.cbegin()) {
						phase_factor( Ntrans, *r, vark, E3 );
						compute_starter( Ntrans, deltak, kzero, kprime, E3, F,
							z, z_source->get(), Zg, *r, phi );
						if (use_turbulence) {
							refract_fluctuations( Ntrans, kzero, *r, z, THETA );
							std::memcpy( GAMMA, THETA, Ntrans * sizeof(double) );
						}
						phi[ 0 ] *= 0.5;
						for (i = 0; i < Ntrans; i++) {
							if (i < Nz) {
								if (use_turbulence) {
									phi[i] *= std::exp(THETA[i] * I);
								}
								TL_mat->set( step_num, i, phi[i] * scaling );
							} else {
								phi[i] = 0.0;
							}
						}

					// subsequent steps
					} else {
						step_num++;
						std::complex<double> cmppow = std::pow(betag/kzero,2.0);
						dr = *r - *(r-1);
						E2 = std::exp( I * dr * kzero
							* (std::sqrt(one-cmppow) - one) );
						phase_factor( Ntrans, dr, vark, E3 );
						NCPA::fft( Ntrans, phi, phitilde );
						NCPA::vector_scale( Ntrans, phitilde, deltaz, phitilde );
						NCPA::circshift( phitilde, Ntrans, -1, phitildeminus );
						NCPA::reverse( phitildeminus, Ntrans, phitildeminus );
						for (i = 0; i < Ntrans; i++) {
							cmppow = kprime[i]/kzero;
							E1 = std::exp( I * dr * kzero
								* (std::sqrt(1.0-std::pow(cmppow,2.0)) - 1.0) );
							E1F = E1 * F[ i ];
							integrand[ i ] = (phitilde[i] + phitildeminus[i])*E1F
								- 2.0 * betag * phitildeminus[i] * (
									(E1F - E2)/(kprime[i]+betag)
								);
						}
						NCPA::ifft( Ntrans, integrand, phi );
						NCPA::vector_multiply( Ntrans, phi, E3, phi );
						NCPA::vector_scale( Ntrans, phi,
							(double)Ntrans * deltak / (2.0*PI), phi );
						if (use_turbulence) {
							refract_fluctuations( Ntrans, kzero, *r, z, THETA );
							NCPA::vector_subtract( Ntrans, THETA, GAMMA, dTHETA );
							std::memcpy( GAMMA, THETA, Ntrans * sizeof(double) );
						}
						phi[ 0 ] *= 0.5;
						for (i = 0; i < Ntrans; i++) {
							if (i < Nz) {
								if (use_turbulence) {
									phi[i] *= std::exp(dTHETA[i] * I);
								}
								TL_mat->set( step_num, i, phi[i] * scaling );
							} else {
								phi[i] = 0.0;
							}
						}
					}
				}

				// how do we save it?
				if (use_suite_mode) {
					NCPA::DenseMatrix<std::complex<double>> *subset
						= new NCPA::DenseMatrix<std::complex<double>>(
							r_indices.size(), z_indices.size() );
					std::vector<double> actual_r, actual_z;
					for (size_t subrow = 0; subrow < r_indices.size(); subrow++) {
						actual_r.push_back( r_vec[ r_indices[ subrow ] ] );
						for (size_t subcol = 0; subcol < z_indices.size(); subcol++) {
							if (subrow == 0) {
								actual_z.push_back( z[ z_indices[ subcol ] ] );
							}
							subset->set( subrow, subcol,
								TL_mat->get(
									r_indices[subrow], z_indices[subcol] ) );
						}
					}
					NCPA::TransmissionLossField *TLF =
						new NCPA::TransmissionLossField( *az, *f,
							actual_r, actual_z, subset );
						TL_.push_back( TLF );
					delete subset;
				} else {
					std::vector<double> z_as_vec( z, z+Nz );
					NCPA::TransmissionLossField *TLF
						= new NCPA::TransmissionLossField( *az, *f,
							r_vec, z_as_vec, TL_mat );
					TL_.push_back( TLF );
				}

				// clean up
				// atm->remove_property("_WC_");
				delete [] E3;
				delete [] phi;
				delete [] THETA;
				delete [] dTHETA;
				delete [] GAMMA;
				delete [] phitilde;
				delete [] phitildeminus;
				delete [] integrand;
				if (use_turbulence) {
					delete turbulence;
				}
				delete TL_mat;

				if (verbose && ((N+1) % 10 == 0)) {
					std::cout << "Completed run " << N+1 << " of "
							  << N_suite << std::endl;
				}
			} // loop over N

			atm->remove_property("_CEFF_");
			atm->remove_property("_ALPHA_");
			delete [] z;
			delete [] vark;
			delete [] F;
			delete [] kprime;

			if (az_vec.size() > 1) {
				output1DTL( tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_1D ), true );
			} else {
				output1DTL( tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_1D ), false );
			}

		} // loop over azimuth
	} // loop over frequency

	return 1;
}

std::string NCPA::GFPESolver::tag_filename( const std::string &base )
		const {
	std::string tagged;
	if (filetag.size() == 0) {
		tagged = base;
	} else {
		tagged = filetag + "." + base;
	}
	return tagged;
}

void NCPA::GFPESolver::compute_starter( size_t Ntrans, double deltak,
	double k0, double *kprime, std::complex<double> *E3, double *F,
	double *z, double z0, std::complex<double> Zhat, double deltar,
	std::complex<double> *&phi ) const {

	size_t i;
	std::complex<double> I( 0.0, 1.0 );
	std::complex<double> *integrand = NCPA::zeros<std::complex<double>>( Ntrans );
	std::complex<double> *result = NCPA::zeros<std::complex<double>>( Ntrans );
	std::complex<double> ks = k0 / Zhat;
	std::complex<double> tmpcmp, factor1, factor2, f2_1, f2_2, R;

	for (i = 0; i < Ntrans; i++) {
		R = (kprime[ i ] - ks) / (kprime[ i ] + ks);
		factor1 = std::exp( -I * z0 * kprime[i] )
			+ R * std::exp( I * kprime[i] * z0 );
		f2_1.real( k0*k0 - kprime[i]*kprime[i] );
		f2_2.real( std::pow(kprime[i]/k0,2.0) );
		factor2 =
			std::exp( I*deltar*std::sqrt(f2_1) )
			/ std::pow(1.0 - f2_2,0.25);
		integrand[ i ] = factor1 * factor2 * F[ i ];
	}

	// take the inverse transform
	NCPA::ifft( Ntrans, integrand, result );

	// assemble return array
	for (i = 0; i < Ntrans; i++) {
		phi[ i ] = std::sqrt(I/(2.0*PI*k0)) * (double)Ntrans * deltak
				* E3[i] * result[i];
		f2_1 = k0*k0 - ks*ks;
		f2_2 = 1.0 - std::pow( ks/k0, 2.0 );
		phi[ i ] -= std::sqrt(2.0*PI/(I*k0)) * 2.0 * ks * std::exp(-I*ks*(z[i]+z0))
				* std::exp(I*std::sqrt(f2_1)*deltar)
				/ std::pow( f2_2, 0.25 );
		phi[ i ] *= std::exp( -I * k0 * deltar );
	}

	delete [] integrand;
	delete [] result;
}


void NCPA::GFPESolver::refract_fluctuations( size_t nz,
		double k_a, double r, double *z, double *&Gamma ) const {

	size_t i, j, nt;
	nt = turbulence->size();

	// build matrices
	// NCPA::Matrix<double> *vec1, *repmat, *sinmat, *kzmat, *mat_Gamma;
	NCPA::Matrix<double> *vec1, *mat1, *mat_Gamma;

	vec1 = new NCPA::DenseMatrix<double>( 1, nt );
	mat1 = new NCPA::DenseMatrix<double>( nt, nz );

	// fill vector and matrix
	for (i = 0; i < nt; i++) {
		vec1->set( 0, i,
			k_a * turbulence->get_G( i ) / turbulence->get_k( i ).real() );
		for (j = 0; j < nz; j++) {
			double temp = r * turbulence->get_k( i ).real()
						+ turbulence->get_alpha( i )
						+ turbulence->get_k( i ).imag() * z[ j ];
			mat1->set( i, j, std::sin( temp ) );
		}
	}

	mat_Gamma = vec1->multiply( mat1 );

	for (j = 0; j < nz; j++) {
		Gamma[ j ] = mat_Gamma->get( 0, j );
	}

	delete vec1;
	delete mat1;
	delete mat_Gamma;
}


void NCPA::GFPESolver::phase_factor( size_t nz, double deltax,
		std::complex<double> *vark, std::complex<double> *&pf ) {
	std::complex<double> I( 0.0, 1.0 );
	for (size_t i = 0; i < nz; i++) {
		pf[ i ] = std::exp( I * deltax * vark[i] );
	}
}


void NCPA::GFPESolver::wavenumber_filter( size_t nz, double *k, double kzero,
	double f1, double f2, double *&F ) const {
	double k1 = f1 * kzero;
	double k2 = f2 * kzero;

	for (size_t i = 0; i < nz; i++) {
		if (std::fabs(k[i]) <= k1) {
			F[ i ] = 1.0;
		} else if (std::fabs(k[i]) > k2) {
			F[ i ] = 0.0;
		} else {
			F[ i ] = 0.5 * (
					1.0 + std::cos(PI * std::fabs(
						(std::fabs(k[i]) - k1)/(k2-k1)
					))
				);
		}
	}
}


double NCPA::GFPESolver::average_value( const std::string &key,
	double range, double low, double high ) const {

	size_t nz = atm->nz( range );
	double *z = NCPA::zeros<double>( nz );
	double *p = NCPA::zeros<double>( nz );
	atm->get_altitude_vector( range, z );
	atm->get_property_vector( range, key, p );
	double avg = 0.0, integral = 0.0, zrange = 0.0;

	// find closest indices
	size_t ind1 = NCPA::find_closest_index( z, nz, low );
	size_t ind2 = NCPA::find_closest_index( z, nz, high );

	if (ind1 == ind2) {
		avg = p[ ind1 ];
	} else {
		integral = NCPA::trapz( ind2 - ind1 + 1, z + ind1, p + ind1 );
		zrange = z[ ind2 ] - z[ ind1 ];
		avg = integral / zrange;
	}
	delete [] z;
	delete [] p;
	return avg;
}

double NCPA::GFPESolver::boundary_layer( double z, double A, double ztop,
		double zatt ) const {
	double alpha = 0.0;
	if (z > zatt) {
		alpha = (A / std::pow( ztop - zatt, 2.0 ))
				* std::pow( z - zatt, 2.0 );
	}
	return alpha;
}

void NCPA::GFPESolver::calculate_effective_sound_speed(
	double az, const std::string &new_key ) {

	// first: was it given explicitly using column "CEFF"?
	if (atm->contains_vector( 0.0, "CEFF" )) {
		atm->convert_property_units( "CEFF", Units::fromString( "m/s" ) );
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

std::complex<double> NCPA::GFPESolver::ground_impedence( double f ) {
	double ratio = sigma / f;
	std::complex<double> Zg(
		1.0 + 0.0511 * std::pow( ratio, 0.75 ),
		0.0768 * std::pow( ratio, 0.73 ) );
	return Zg;
}


void NCPA::GFPESolver::output1DTL( double height,
		const std::string &filename, bool append ) const {

	if (f_vec.size() > 1) {
		throw std::runtime_error( "Undefined behavior: no format for output1DTL for multiple frequencies" );
	}

	std::ofstream out_1d;
	if (append) {
		out_1d.open( filename, std::ofstream::out | std::ofstream::app );
		out_1d << std::endl;
	} else {
		out_1d.open( filename, std::ofstream::out | std::ofstream::trunc );
	}

	// iterate through transmission loss fields.  For now, assume they're
	// in azimuth order (they should be)
	std::vector<NCPA::TransmissionLossField *>::const_iterator it;
	for (it = TL_.cbegin(); it != TL_.cend(); ++it) {
		double az, freq;
		// size_t nr, nz;
		// double *r, *z;
		// std::complex<double> **tl;
		// (*it)->as_arrays( az, freq, nr, r, nz, z, tl );
		std::vector<double> r, z;
		NCPA::DenseMatrix<std::complex<double>> *tl;
		(*it)->as_vectors( az, freq, r, z, tl );

		// find height index
		// size_t zi = NCPA::find_closest_index<double>( z, nz, height );
		size_t zi = NCPA::find_closest_index<double>( z, height );
		for (size_t i = 0; i < r.size(); i++) {
			out_1d << r[ i ]/1000.0 << " " << az << " "
				   << tl->get( i, zi ).real() << " "
				   << tl->get( i, zi ).imag() << std::endl;
		}

		// more?
		if (!use_suite_mode && ((it + 1) != TL_.cend())) {
			out_1d << std::endl;
		}

		delete tl;
	}
	out_1d.close();
}

void NCPA::GFPESolver::get1DTL( size_t n, double height,
		std::vector<double> &r,
		std::vector<std::complex<double>> &tl_vec ) const {
	std::vector<double> z;
	double tl_az, tl_freq;
	NCPA::DenseMatrix<std::complex<double>> *tl;
	r.clear();
	tl_vec.clear();
	TL_[n]->as_vectors( tl_az, tl_freq, r, z, tl );

	size_t zi = NCPA::find_closest_index<double>( z, height );
	tl->get_column( zi, tl_vec );
	delete tl;
}

void NCPA::GFPESolver::get1DTL( double height, std::vector<double> &r,
		std::vector<std::complex<double>> &tl_vec ) const {
	size_t n = 0;
	this->get1DTL( n, height, r, tl_vec );
}

void NCPA::GFPESolver::get1DTL( size_t n, std::vector<double> &r,
		std::vector<std::complex<double>> &tl_vec ) const {
	this->get1DTL( n, z_receiver->get(), r, tl_vec );
}

void NCPA::GFPESolver::get1DTL( std::vector<double> &r,
		std::vector<std::complex<double>> &tl_vec ) const {
	size_t n = 0;
	this->get1DTL( n, z_receiver->get(), r, tl_vec );
}

void NCPA::GFPESolver::output1DTL( bool append ) const {
	std::string fn = tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_1D );
	output1DTL( fn, append );
}

void NCPA::GFPESolver::output1DTL( double height, bool append ) const {
	std::string fn = tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_1D );
	output1DTL( height, fn, append );
}

void NCPA::GFPESolver::output1DTL( const std::string &filename,
		bool append ) const {
	output1DTL( z_receiver->get(), filename, append );
}

void NCPA::GFPESolver::output2DTL( bool append ) const {
	std::string fn = tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_2D );
	output2DTL( fn, append );
}

void NCPA::GFPESolver::output2DTL( const std::string &filename,
			bool append ) const {
	this->output2DTL( 0, filename, append );
}

void NCPA::GFPESolver::output2DTL( size_t n, bool append ) const {
	std::string filename = tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_2D );
	this->output2DTL( n, filename, append );
}

void NCPA::GFPESolver::output2DTL( size_t n, const std::string &filename,
			bool append ) const {

	if (f_vec.size() > 1) {
		throw std::runtime_error( "Undefined behavior: no format for output2DTL for multiple frequencies" );
	}

	std::ofstream out_2d;
	if (append) {
		out_2d.open( filename, std::ofstream::out | std::ofstream::app );
	} else {
		out_2d.open( filename, std::ofstream::out | std::ofstream::trunc );
	}
	double az, freq;
	std::vector<double> r, z;
	NCPA::DenseMatrix<std::complex<double>> *tl;
	TL_[n]->as_vectors( az, freq, r, z, tl );

	for (size_t i = 0; i < r.size(); i++) {
		for (size_t j = 0; j < z.size(); j++) {
			out_2d << az << " "
				   << r[ i ]/1000.0 << " "
				   << z[ j ]/1000.0 << " "
				   << tl->get( i, j ).real() << " "
				   << tl->get( i, j ).imag() << std::endl;
		}
		if (!use_suite_mode) {
			out_2d << std::endl;
		}
	}
	out_2d.close();
	delete tl;
}

void NCPA::GFPESolver::get2DTL( size_t n, std::vector<double> &r,
		std::vector<double> &z,
		NCPA::DenseMatrix<std::complex<double>> *&tl_mat ) const {

	double az, freq;
	TL_[n]->as_vectors( az, freq, r, z, tl_mat );
}

void NCPA::GFPESolver::get2DTL( std::vector<double> &r,
		std::vector<double> &z,
		NCPA::DenseMatrix<std::complex<double>> *&tl_mat ) const {

	double az, freq;
	TL_[0]->as_vectors( az, freq, r, z, tl_mat );
}

size_t NCPA::GFPESolver::runs() const {
	return TL_.size();
}

void NCPA::GFPESolver::truncateFile( const std::string &filename ) const {
	std::ofstream tfile( filename,
		std::ofstream::out | std::ofstream::trunc);
	tfile.close();
}

void NCPA::GFPESolver::truncate1DFile() const {
	std::string filename = tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_1D );
	truncateFile( filename );
}

void NCPA::GFPESolver::truncate2DFile() const {
	std::string filename = tag_filename( NCPAPROP_GFPE_TLOSS_FILENAME_2D );
	truncateFile( filename );
}
