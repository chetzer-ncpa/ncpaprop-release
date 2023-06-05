#include "Atmosphere2D.h"
#include "Atmosphere1D.h"
#include <vector>
#include <climits>
#include <cfloat>
#include <stdexcept>
#include <algorithm>
#include <cstring>

#include "gsl/gsl_version.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"



NCPA::Atmosphere2D::Atmosphere2D() {
	profiles_.clear();
	midpoints_.clear();
	clear_last_index_();
	range_units_ = NCPA::UNITS_NONE;
	max_valid_range_ = 0;
	
	topo_ground_heights_ = NULL;
	topo_ranges_ = NULL;
	topo_accel_ = NULL;
	topo_spline_ = NULL;
}

NCPA::Atmosphere2D::~Atmosphere2D() {
	free_ground_elevation_spline_();
	for ( std::vector< NCPA::Atmosphere1D * >::iterator i = profiles_.begin();
		i != profiles_.end(); ++i ) {
		delete (*i);
	}
	profiles_.clear();
	midpoints_.clear();
}

void NCPA::Atmosphere2D::set_maximum_valid_range( double maxrange ) {
	if (range_units_ == NCPA::UNITS_NONE) {
		throw std::runtime_error( "Range units have not been specified!" );
	}

	max_valid_range_ = maxrange;

	std::cout << "Setting maximum range to " << maxrange << " " << NCPA::Units::toString( range_units_ ) << std::endl;
}

double NCPA::Atmosphere2D::get_maximum_valid_range() const {
	return max_valid_range_;
}

bool NCPA::sort_profiles_by_range_( Atmosphere1D *p1, Atmosphere1D *p2 ) {
	// assumes consistent units
	return p1->get( "_RANGE_" ) < p2->get( "_RANGE_" );
}

void NCPA::Atmosphere2D::insert_profile( const NCPA::Atmosphere1D *profile, double range ) {

	// throw an exception if units haven't been specified yet
	if ( range_units_ == NCPA::UNITS_NONE ) {
		throw std::runtime_error( "Range units have not been specified!" );
	}
	NCPA::Atmosphere1D *newProfile = new NCPA::Atmosphere1D( *profile );
	newProfile->add_property( "_RANGE_", range, range_units_ );
	profiles_.push_back( newProfile );
	sorted_ = false;
	clear_last_index_();

	if (range > max_valid_range_) {
		set_maximum_valid_range( range );
	}

	// clear ground elevation spline
	free_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::set_insert_range_units( units_t u ) {
	range_units_ = u;
	free_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::clear_last_index_() {
	last_index_ = INT_MAX;
	last_index_min_range_ = 0.0;
	last_index_max_range_ = 0.0;
}

void NCPA::Atmosphere2D::set_last_index_( size_t ind ) {
	if (profiles_.size() == 0) {
		// should never happen, but you never know
		throw std::runtime_error( "No profiles have been added to 2-D atmosphere!" );
	} else if (profiles_.size() == 1) {
		last_index_ = 0;
		last_index_min_range_ = 0.0;
		last_index_max_range_ = DBL_MAX;
	} else if (ind == 0) {
		last_index_ = 0;
		last_index_min_range_ = 0.0;
		last_index_max_range_ = midpoints_[ 0 ];
	} else if (ind == (profiles_.size() - 1) ) {
		last_index_ = ind;
		last_index_min_range_ = midpoints_[ ind-1 ];
		last_index_max_range_ = DBL_MAX;
	} else {
		last_index_ = ind;
		last_index_min_range_ = midpoints_[ ind-1 ];
		last_index_max_range_ = midpoints_[ ind ];
	}
}

void NCPA::Atmosphere2D::sort_profiles() {
	std::sort( profiles_.begin(), profiles_.end(), NCPA::sort_profiles_by_range_ );
	midpoints_.clear();
	clear_last_index_();
	sorted_ = true;
	free_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::calculate_midpoints_() {
	midpoints_.clear();
	if (profiles_.size() == 1) {
		return;
	}

	if (!sorted_) {
		sort_profiles();
	}

	for (size_t i = 1; i < profiles_.size(); i++) {
		double r0 = profiles_.at( i-1 )->get( "_RANGE_" );
		double r1 = profiles_.at( i )->get( "_RANGE_" );
		double mid = r0 + ( r1 - r0 ) * 0.5;
		midpoints_.push_back( mid );
	}
}

size_t NCPA::Atmosphere2D::get_profile_index( double range ) {
	if (profiles_.size() == 0) {
		throw std::runtime_error( "No profiles have been added to 2-D atmosphere!" );
	}

	if (profiles_.size() == 1) {
		return 0;
	}

	if (midpoints_.size() == 0) {
		calculate_midpoints_();
		clear_last_index_();
	}

	// first check to see if we have a valid last range, to avoid going through the search again
	if (last_index_ < INT_MAX) {
		if (range >= last_index_min_range_ && range < last_index_max_range_) {
			return last_index_;
		}
	}

	size_t ind = 0;
	while ( ind < midpoints_.size() && midpoints_[ ind ] <= range ) {
		ind++;
	}
	set_last_index_( ind );
	return ind;
}

double NCPA::Atmosphere2D::get( double range, const std::string &key ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get( key );
}

double NCPA::Atmosphere2D::get( double range, const std::string &key, double altitude ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get( key, altitude );
}

double NCPA::Atmosphere2D::get_first_derivative( double range, const std::string &key,
		double altitude ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get_first_derivative( key, altitude );
}

double NCPA::Atmosphere2D::get_second_derivative( double range, const std::string &key,
		double altitude ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get_second_derivative( key, altitude );
}

size_t NCPA::Atmosphere2D::nz( double range ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->nz();
}

void NCPA::Atmosphere2D::get_altitude_vector( double range, double *buffer,
		units_t *buffer_units ) {
	size_t ind = get_profile_index( range );
	profiles_.at( ind )->get_altitude_vector( buffer, buffer_units );
}

void NCPA::Atmosphere2D::get_property_vector( double range, const std::string &key,
	double *buffer, units_t *buffer_units ) {
	size_t ind = get_profile_index( range );
	profiles_.at( ind )->get_property_vector( key, buffer, buffer_units );
}

void NCPA::Atmosphere2D::get_altitude_vector( double range, double *buffer ) {
	size_t ind = get_profile_index( range );
	profiles_.at( ind )->get_altitude_vector( buffer );
}

void NCPA::Atmosphere2D::get_property_vector( double range, const std::string &key,
		double *buffer ) {
	size_t ind = get_profile_index( range );
	profiles_.at( ind )->get_property_vector( key, buffer );
}

NCPA::units_t NCPA::Atmosphere2D::get_range_units() const {
	return range_units_;
}

NCPA::units_t NCPA::Atmosphere2D::get_altitude_units( double range ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get_altitude_units();
}

NCPA::units_t NCPA::Atmosphere2D::get_property_units( double range,
		const std::string &key ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get_property_units( key );
}

double NCPA::Atmosphere2D::get_minimum_altitude( double range ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get_minimum_altitude();
}

double NCPA::Atmosphere2D::get_maximum_altitude( double range ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->get_maximum_altitude();
}

bool NCPA::Atmosphere2D::contains_scalar( double range, const std::string &key ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->contains_scalar( key );
}

bool NCPA::Atmosphere2D::contains_vector( double range, const std::string &key ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->contains_vector( key );
}

bool NCPA::Atmosphere2D::contains_key( double range, const std::string &key ) {
	size_t ind = get_profile_index( range );
	return profiles_.at( ind )->contains_key( key );
}

void NCPA::Atmosphere2D::calculate_density_from_temperature_and_pressure(
		const std::string &new_key, const std::string &temperature_key,
		const std::string &pressure_key, units_t density_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_density_from_temperature_and_pressure(
			new_key, temperature_key, pressure_key, density_units );
	}
}

void NCPA::Atmosphere2D::calculate_sound_speed_from_temperature( const std::string &new_key, const std::string &temperature_key,
		units_t wind_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_sound_speed_from_temperature( new_key, temperature_key, wind_units );
	}
}

void NCPA::Atmosphere2D::calculate_sound_speed_from_pressure_and_density( const std::string &new_key, const std::string &pressure_key,
		const std::string &density_key, NCPA::units_t wind_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_sound_speed_from_pressure_and_density( new_key, pressure_key, density_key, wind_units );
	}
}

void NCPA::Atmosphere2D::calculate_wind_speed( const std::string &new_key, const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_wind_speed( new_key, we_wind_speed_key, sn_wind_speed_key );
	}
}

void NCPA::Atmosphere2D::calculate_wind_direction( const std::string &new_key, const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
		NCPA::units_t direction_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_wind_direction( new_key, we_wind_speed_key, sn_wind_speed_key, direction_units );
	}
}

void NCPA::Atmosphere2D::calculate_attenuation( const std::string &new_key,
		const std::string &temperature_key, const std::string &pressure_key,
		const std::string &density_key, double freq, double tweak_factor ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_attenuation( new_key, temperature_key,
			pressure_key, density_key, freq, tweak_factor );
	}
}

void NCPA::Atmosphere2D::calculate_attenuation( const std::string &new_key,
		const std::string &temperature_key, const std::string &pressure_key,
		const std::string &density_key, const std::string &humidity_key,
		double freq, double tweak_factor ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_attenuation( new_key, temperature_key,
			pressure_key, density_key, humidity_key, freq, tweak_factor );
	}
}


void NCPA::Atmosphere2D::read_attenuation_from_file( const std::string &new_key, const std::string &filename ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->read_attenuation_from_file( new_key, filename );
	}
}

void NCPA::Atmosphere2D::calculate_wind_component( const std::string &new_key, const std::string &wind_speed_key, const std::string &wind_direction_key,
		double azimuth ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_wind_component( new_key, wind_speed_key, wind_direction_key, azimuth );
	}
}

void NCPA::Atmosphere2D::calculate_effective_sound_speed( const std::string &new_key, const std::string &sound_speed_key,
	const std::string &wind_component_key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_effective_sound_speed( new_key, sound_speed_key, wind_component_key );
	}
}

void NCPA::Atmosphere2D::convert_altitude_units( NCPA::units_t new_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->convert_altitude_units( new_units );
	}
	free_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::convert_property_units( const std::string &key, NCPA::units_t new_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->convert_property_units( key, new_units );
	}
	if (key == "Z0") {
		free_ground_elevation_spline_();
	}
}

void NCPA::Atmosphere2D::convert_range_units( NCPA::units_t new_units ) {
	NCPA::units_t old_units = range_units_;
	this->convert_property_units( "_RANGE_", new_units );
	max_valid_range_ = NCPA::Units::convert( max_valid_range_, old_units, new_units );
	range_units_ = new_units;
	sorted_ = false;
	clear_last_index_();
	calculate_midpoints_();
	free_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::get_minimum_altitude_limits( double &lowlimit, double &highlimit ) {
	lowlimit = DBL_MAX;
	highlimit = -DBL_MAX;
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		double curalt = (*it)->get_minimum_altitude();
		highlimit = NCPA::max( highlimit, curalt );
		lowlimit = NCPA::min( lowlimit, curalt );
	}
}

void NCPA::Atmosphere2D::get_maximum_altitude_limits( double &lowlimit, double &highlimit ) {
	lowlimit = DBL_MAX;
	highlimit = -DBL_MAX;
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		double curalt = (*it)->get_maximum_altitude();
		highlimit = NCPA::max( highlimit, curalt );
		lowlimit = NCPA::min( lowlimit, curalt );
	}
}


void NCPA::Atmosphere2D::add_property( const std::string &key, size_t n_points, double *quantity_points,
			units_t quantity_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->add_property( key, n_points, quantity_points, quantity_units );
	}
}

void NCPA::Atmosphere2D::add_property( const std::string &key, double value, units_t units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->add_property( key, value, units );
	}
}

void NCPA::Atmosphere2D::copy_vector_property( const std::string &old_key, const std::string &new_key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->copy_vector_property( old_key, new_key );
	}
}

void NCPA::Atmosphere2D::copy_scalar_property( const std::string &old_key, const std::string &new_key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->copy_scalar_property( old_key, new_key );
	}
	if (old_key == "Z0") {
		free_ground_elevation_spline_();
	}
}

void NCPA::Atmosphere2D::remove_property( const std::string &key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->remove_property( key );
	}
	if (key == "Z0") {
		free_ground_elevation_spline_();
	}
}

std::vector< NCPA::Atmosphere1D * >::iterator NCPA::Atmosphere2D::first_profile() {
	return profiles_.begin();
}

std::vector< NCPA::Atmosphere1D * >::iterator NCPA::Atmosphere2D::last_profile() {
	return profiles_.end();
}

void NCPA::Atmosphere2D::read_elevation_from_file( const std::string &filename ) {
	override_profile_z0_ = true;
	ground_elevation_file_ = filename;
	generate_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::finalize_elevation_from_profiles() {
	override_profile_z0_ = false;
	generate_ground_elevation_spline_();
}

void NCPA::Atmosphere2D::generate_ground_elevation_spline_() {

	free_ground_elevation_spline_();

	topo_accel_ = gsl_interp_accel_alloc();
	if (override_profile_z0_) {
		setup_ground_elevation_spline_from_file_();
	} else {
		setup_ground_elevation_spline_from_profiles_();
	}
}

void NCPA::Atmosphere2D::setup_ground_elevation_spline_from_file_() {

	std::cout << "Reading topography from " << ground_elevation_file_
			  << std::endl;
	std::ifstream topofile( ground_elevation_file_ );
	units_t z_units = get_altitude_units( 0.0 );
	units_t r_units = range_units_;
	units_t file_r_units = NCPA::Units::fromString( "km" );
	units_t file_z_units = NCPA::Units::fromString( "m" );
	std::string line;
	std::string delims = ":,= ";
	std::vector< double > rvec, zvec;

	std::getline( topofile, line );
	while( topofile.good() ) {
		line = NCPA::deblank( line );
		if (line.size() > 0) {
			if (line[ 0 ] == '#') {
				if (line.size() > 1 && line[ 1 ] == '%') {

					line.erase(0,2);
					line = NCPA::deblank( line );
					std::vector<std::string> fields = NCPA::split(
						line, delims );
					if (fields.size() < 2) {
						std::cerr << "Topographic file descriptive header line "
								  << line << " has no delimiter characters ("
								  << delims << "), ignoring" << std::endl;
					} else {
						units_t tempunits = NCPA::Units::fromString( fields[1] );
						if (tempunits == UNITS_NONE) {
							std::cerr << "Unrecognized units " << fields[1]
									  << ", ignoring" << std::endl;
						} else {
							switch ((fields[0])[0]) {
								case 'r':
								case 'R':
									file_r_units = tempunits;
									break;
								case 'z':
								case 'Z':
									file_z_units = tempunits;
									break;
								default:
									std::cerr << "Unrecognized parameter tag " << line[0]
											  << ", must be in [RrZz].  Ignoring" << std::endl;
							}
						}
					}


					// size_t delimpos = line.find_last_of( delims );
					// if (delimpos == std::string::npos) {
					// 	std::cerr << "Topographic file descriptive header line "
					// 			  << line << " has no delimiter characters ("
					// 			  << delims << "), ignoring" << std::endl;
					// } else {
					// 	// chop off first two characters
					// 	line.erase(0,2);
					// 	line = NCPA::deblank( line );
					// 	delims += " ";
					// 	delimpos = line.find_last_of( delims );
					// 	std::string ustr = NCPA::deblank(line.substr( delimpos+1 ));
					// 	units_t tempunits = NCPA::Units::fromString( ustr );
					// 	if (tempunits == UNITS_NONE) {
					// 		std::cerr << "Unrecognized units " << ustr << ", ignoring"
					// 				  << std::endl;
					// 	} else {
					// 		switch (line[0]) {
					// 			case 'r':
					// 			case 'R':
					// 				file_r_units = tempunits;
					// 				break;
					// 			case 'z':
					// 			case 'Z':
					// 				file_z_units = tempunits;
					// 				break;
					// 			default:
					// 				std::cerr << "Unrecognized parameter tag " << line[0]
					// 						  << ", must be in [RrZz].  Ignoring" << std::endl;
					// 		}
					// 	}
					// }
				}
			} else {
				std::vector< std::string > parts = NCPA::split( line );
				double r = std::stod( parts[ 0 ] );
				double z = std::stod( parts[ 1 ] );
				rvec.push_back( NCPA::Units::convert( r, file_r_units, r_units ) );
				zvec.push_back( NCPA::Units::convert( z, file_z_units, z_units ) );
			}
		}

		std::getline( topofile, line );
	}
	topofile.close();

	// Now create the splines
	size_t np = rvec.size();
#if GSL_MAJOR_VERSION > 1
	topo_spline_ = gsl_spline_alloc( gsl_interp_steffen, np + 2 );
#else
	std::cout << "Version 1 of GSL does not include Steffen interpolation."
			  << std::endl
			  << "Defaulting to Cubic Spline interpolation for ground height."
			  << std::endl;
	topo_spline_ = gsl_spline_alloc( gsl_interp_cspline, np + 2 );
#endif
	topo_ground_heights_ = new double[ np + 2 ];
	topo_ranges_ = new double[ np + 2 ];
	std::memset( topo_ground_heights_, 0, (np + 2) * sizeof(double) );
	std::memset( topo_ranges_, 0, (np + 2) * sizeof(double) );

	topo_ground_heights_[ 0 ] = rvec.front();
	topo_ranges_[ 0 ] = NCPA::Units::convert( -1000.0, UNITS_DISTANCE_METERS, r_units );
	for (size_t i = 0; i < np; i++) {
		topo_ground_heights_[ i+1 ] = zvec[ i ];
		topo_ranges_[ i+1 ] = rvec[ i ];
	}
	topo_ground_heights_[ np+1 ] = topo_ground_heights_[ np ];
	topo_ranges_[ np+1 ] = NCPA::max<double>(
		topo_ranges_[ np ]
			+ NCPA::Units::convert( 1000.0, UNITS_DISTANCE_METERS, r_units ),
		NCPA::Units::convert( get_maximum_valid_range(),
			get_range_units(), r_units )
		);

	gsl_spline_init( topo_spline_, topo_ranges_, topo_ground_heights_, np+2 ); 
}

void NCPA::Atmosphere2D::setup_ground_elevation_spline_from_profiles_() {
	size_t np = profiles_.size();
#if GSL_MAJOR_VERSION > 1
	topo_spline_ = gsl_spline_alloc( gsl_interp_steffen, np + 2 );
#else
	std::cout << "Version 1 of GSL does not include Steffen interpolation."
			  << std::endl
			  << "Defaulting to Cubic Spline interpolation for ground height."
			  << std::endl;
	topo_spline_ = gsl_spline_alloc( gsl_interp_cspline, np + 2 );
#endif
	topo_ground_heights_ = new double[ np + 2 ];
	topo_ranges_ = new double[ np + 2 ];
	std::memset( topo_ground_heights_, 0, (np + 2) * sizeof(double) );
	std::memset( topo_ranges_, 0, (np + 2) * sizeof(double) );

	std::vector< NCPA::Atmosphere1D * >::iterator it = this->first_profile();
	topo_ground_heights_[ 0 ] = (*it)->get( "Z0" );
	topo_ranges_[ 0 ] = NCPA::Units::convert( -1000.0, UNITS_DISTANCE_METERS, range_units_ );
	size_t pnum = 1;
	for ( ; it != this->last_profile(); ++it ) {
		if ((*it)->contains_scalar("Z0")) {
			topo_ground_heights_[ pnum ] = (*it)->get( "Z0" );
		} else {
			topo_ground_heights_[ pnum ] = 0.0;
		}
		topo_ranges_[ pnum++ ] = (*it)->get( "_RANGE_" );
	}
	topo_ground_heights_[ pnum ] = topo_ground_heights_[ pnum-1 ];
	topo_ranges_[ pnum ] = NCPA::Units::convert( 1000.0, UNITS_DISTANCE_METERS, range_units_ )
		+ max_valid_range_;

	gsl_spline_init( topo_spline_, topo_ranges_, topo_ground_heights_, np+2 );
}

void NCPA::Atmosphere2D::free_ground_elevation_spline_() {
	if (topo_spline_ != NULL) {
		gsl_spline_free( topo_spline_ );
		topo_spline_ = NULL;
	}
	if (topo_accel_ != NULL) {
		gsl_interp_accel_free( topo_accel_ );
		topo_accel_ = NULL;
	}
	if (topo_ranges_ != NULL) {
		delete [] topo_ranges_;
		topo_ranges_ = NULL;
	}
	if (topo_ground_heights_ != NULL) {
		delete [] topo_ground_heights_;
		topo_ground_heights_ = NULL;
	}
}

double NCPA::Atmosphere2D::get_interpolated_ground_elevation( double range ) {
	if (topo_spline_ == NULL) {
		generate_ground_elevation_spline_();
	}
	return gsl_spline_eval( topo_spline_, range, topo_accel_ );
}

double NCPA::Atmosphere2D::get_interpolated_ground_elevation_first_derivative( double range ) {
	if (topo_spline_ == NULL) {
		generate_ground_elevation_spline_();
	}
	return gsl_spline_eval_deriv( topo_spline_, range, topo_accel_ );
}

double NCPA::Atmosphere2D::get_interpolated_ground_elevation_second_derivative( double range ) {
	if (topo_spline_ == NULL) {
		generate_ground_elevation_spline_();
	}
	return gsl_spline_eval_deriv2( topo_spline_, range, topo_accel_ );
}


void NCPA::Atmosphere2D::print_atmosphere(
			const std::vector< std::string >& columnorder,
			double range, const std::string &altitude_key,
			std::ostream& os ) {

	// check columnorder variable for key validity
	std::vector< std::string >::const_iterator vit;
	for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit) {
		if (! contains_vector( range, *vit ) ) {
			throw std::invalid_argument( "No vector quantity exists with key " + *vit );
		}
	}

	// Now column descriptors.  Altitude first
	NCPA::units_t z_units;
	size_t nz_ = nz( range );
	double *z_ = NCPA::zeros<double>( nz_ );
	get_altitude_vector( range, z_, &z_units );
	os  << "#% 1, " << altitude_key << ", "
		<< NCPA::Units::toStr( z_units ) << std::endl;
	unsigned int column = 2;
	for ( vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
		os  << "#% " << column << ", "
			<< remove_underscores( *vit ) << ", "
			<< NCPA::Units::toStr( get_property_units( range, *vit ) )
			<< std::endl;
		column++;
	}

	// Now columns
	os.setf( std::ios::scientific, 	std::ios::floatfield );
	os.setf( std::ios::right, 		std::ios::adjustfield );
	os.precision( 6 );
	os.width( 9 );
	os.fill( ' ' );
	for ( size_t i = 0; i < nz_; i++) {
		os << z_[ i ];
		for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
			os << " " << get( range, *vit, z_[ i ] );
		}
		os << std::endl;
	}
	os.flush();

	delete [] z_;
}