#include "NCPACommon.h"
#include "StratifiedAtmosphere3D.h"

#include <iostream>
#include <cfloat>
#include <stdexcept>


NCPA::StratifiedAtmosphere3D::StratifiedAtmosphere3D( const Atmosphere1D *atm ) {
	profile_ = new NCPA::Atmosphere1D( *atm );
	range_units_ = UNITS_DISTANCE_KILOMETERS;
}

NCPA::StratifiedAtmosphere3D::StratifiedAtmosphere3D( const std::string &filename,
	std::string headerfilename ) {
	profile_ = new NCPA::Atmosphere1D( filename, headerfilename );
	range_units_ = UNITS_DISTANCE_KILOMETERS;
}

NCPA::StratifiedAtmosphere3D::~StratifiedAtmosphere3D() {
	delete profile_;
}

void NCPA::StratifiedAtmosphere3D::convert_range_units( units_t new_units ) {
	range_units_ = new_units;
}

void NCPA::StratifiedAtmosphere3D::convert_altitude_units( units_t new_units ) {
	profile_->convert_altitude_units( new_units );
}

void NCPA::StratifiedAtmosphere3D::convert_property_units( const std::string &key,
	units_t new_units ) {
	profile_->convert_property_units( key, new_units );
}

NCPA::units_t NCPA::StratifiedAtmosphere3D::get_range_units() const {
	return range_units_;
}

NCPA::units_t NCPA::StratifiedAtmosphere3D::get_altitude_units() const {
	return profile_->get_altitude_units();
}

NCPA::units_t NCPA::StratifiedAtmosphere3D::get_property_units( const std::string &key ) const {
	return profile_->get_property_units( key );
}


// querying
double NCPA::StratifiedAtmosphere3D::get_maximum_valid_range( double azimuth ) const {
	return DBL_MAX;
}

bool NCPA::StratifiedAtmosphere3D::contains_property( const std::string &key ) const {
	return profile_->contains_key( key );
}

bool NCPA::StratifiedAtmosphere3D::contains_vector( const std::string &key ) const {
	return profile_->contains_vector( key );
}

bool NCPA::StratifiedAtmosphere3D::contains_scalar( const std::string &key ) const {
	return profile_->contains_scalar( key );
}

// data retrieval
double NCPA::StratifiedAtmosphere3D::get( double x, double y, const std::string &key ) const {
	return profile_->get( key );
}

double NCPA::StratifiedAtmosphere3D::get( double x, double y, double z,
			const std::string &key ) const {
	return profile_->get( key, z );
}

double NCPA::StratifiedAtmosphere3D::get_derivative( double x, double y,
	const std::string &key, size_t order, deriv_t *directions ) const {
	return 0.0;
}

double NCPA::StratifiedAtmosphere3D::get_derivative( double x, double y,
			const std::string &key, deriv_t direction ) const {
	return 0.0;
}

double NCPA::StratifiedAtmosphere3D::get_derivative( double x, double y, double z,
	const std::string &key, size_t order, deriv_t *directions ) const {

	// if requesting only vertical derivatives, return from the profile, otherwise zero
	bool all_z = true;
	for (size_t i = 0; i < order; i++) {
		all_z = all_z && (directions[i] == 2 );
	}
	if (all_z) {
		switch (order) {
			case 1:
				return profile_->get_first_derivative( key, z );
				break;
			case 2:
				return profile_->get_second_derivative( key, z );
				break;
			default:
				throw std::runtime_error( "Maximum second derivative available" );
		}
	}

	return 0.0;
}

void NCPA::StratifiedAtmosphere3D::add_property( const std::string &key, double ***prop,
			size_t nx, size_t ny, size_t nz, NCPA::units_t units ) {
	profile_->add_property( key, nz, prop[ 0 ][ 0 ], units );
}

void NCPA::StratifiedAtmosphere3D::add_property( const std::string &key, double **prop,
			size_t nx, size_t ny, NCPA::units_t units ) {
	profile_->add_property( key, prop[ 0 ][ 0 ], units );
}

void NCPA::StratifiedAtmosphere3D::calculate_sound_speed_from_temperature(
			const std::string &new_key, const std::string &temperature_key,
			NCPA::units_t wind_units ) {
	profile_->calculate_sound_speed_from_temperature( new_key, temperature_key,
		wind_units );
}

void NCPA::StratifiedAtmosphere3D::calculate_sound_speed_from_pressure_and_density(
			const std::string &new_key, const std::string &pressure_key,
			const std::string &density_key, units_t wind_units ) {
	profile_->calculate_sound_speed_from_pressure_and_density( new_key,
		pressure_key, density_key, wind_units );
}

void NCPA::StratifiedAtmosphere3D::calculate_wind_speed( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key ) {
	profile_->calculate_wind_speed( new_key, we_wind_speed_key, sn_wind_speed_key );
}

void NCPA::StratifiedAtmosphere3D::calculate_wind_direction( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
			units_t direction_units ) {
	profile_->calculate_wind_direction( new_key, we_wind_speed_key, sn_wind_speed_key,
		direction_units );
}

void NCPA::StratifiedAtmosphere3D::calculate_attenuation( const std::string &new_key,
			const std::string &temperature_key, const std::string &pressure_key,
			const std::string &density_key, double freq, double tweak_factor ) {
	profile_->calculate_attenuation( new_key, temperature_key, pressure_key,
		density_key, freq, tweak_factor );
}

void NCPA::StratifiedAtmosphere3D::calculate_wind_component( const std::string &new_key,
			const std::string &wind_speed_key, const std::string &wind_direction_key,
			double azimuth ) {
	profile_->calculate_wind_component( new_key, wind_speed_key, wind_direction_key,
		azimuth );
}

void NCPA::StratifiedAtmosphere3D::calculate_effective_sound_speed(
	const std::string &new_key, const std::string &sound_speed_key,
	const std::string &wind_component_key ) {
	profile_->calculate_effective_sound_speed( new_key, sound_speed_key, wind_component_key );
}


void NCPA::StratifiedAtmosphere3D::remove_property( const std::string &key ) {
	profile_->remove_property( key );
}

void NCPA::StratifiedAtmosphere3D::remove_vector_property( const std::string &key ) {
	profile_->remove_property( key );
}

void NCPA::StratifiedAtmosphere3D::remove_scalar_property( const std::string &key ) {
	profile_->remove_property( key );
}

double NCPA::StratifiedAtmosphere3D::get_minimum_altitude( double x, double y ) const {
	return profile_->get_minimum_altitude();
}

void NCPA::StratifiedAtmosphere3D::get_property_template( const std::string &basis,
				size_t &nx, double *x, NCPA::units_t &x_units,
				size_t &ny, double *y, NCPA::units_t &y_units,
				double **&prop ) const {

	prop = NCPA::allocate_matrix<double>( 1, 1 );
	nx = 1;
	x = new double[ 1 ];
	x[ 0 ] = 0.0;
	x_units = NCPA::UNITS_DISTANCE_KILOMETERS;
	ny = 1;
	y = new double[ 1 ];
	y[ 0 ] = 0.0;
	y_units = NCPA::UNITS_DISTANCE_KILOMETERS;



}

void NCPA::StratifiedAtmosphere3D::get_property_template( const std::string &basis,
			size_t &nx, double *x, NCPA::units_t &x_units,
			size_t &ny, double *y, NCPA::units_t &y_units,
			size_t &nz, double *z, NCPA::units_t &z_units,
			double ***&prop ) const {
	if (!profile_->contains_vector( basis )) {
		throw std::runtime_error( "Key " + basis + " not found");
	}

	nx = 1;
	x = new double[ 1 ];
	x[ 0 ] = 0.0;
	x_units = NCPA::UNITS_DISTANCE_KILOMETERS;
	ny = 1;
	y = new double[ 1 ];
	y[ 0 ] = 0.0;
	y_units = NCPA::UNITS_DISTANCE_KILOMETERS;
	nz = profile_->nz();
	z = new double[ nz ];
	profile_->get_altitude_vector( z, &z_units );
	prop = NCPA::matrix3d<double>( nx, ny, nz );

}

void NCPA::StratifiedAtmosphere3D::free_property_template(
				size_t nx, double *x,
				size_t ny, double *y,
				double **prop ) const {
	NCPA::free_matrix<double>( prop, nx, ny );
	delete [] x;
	delete [] y;
}

void NCPA::StratifiedAtmosphere3D::free_property_template(
				size_t nx, double *x,
				size_t ny, double *y,
				size_t nz, double *z,
				double ***prop ) const {
	NCPA::free_matrix3d<double>( prop, nx, ny, nz );
	delete [] x;
	delete [] y;
	delete [] z;
}

void NCPA::StratifiedAtmosphere3D::copy_property( const std::string &old_key,
				const std::string &new_key ) {
	if (contains_property(new_key)) {
		throw std::runtime_error( "Requested new key " + new_key + " already exists" );
	}
	if (contains_vector(old_key)) {
		profile_->copy_vector_property( old_key, new_key );
	} else if (contains_scalar(old_key)) {
		profile_->copy_scalar_property( old_key, new_key );
	}
	throw std::runtime_error( "Key " + old_key + " not found" );
}

void NCPA::StratifiedAtmosphere3D::read_attenuation_from_file( const std::string &key,
			const std::string &filename ) {
	profile_->read_attenuation_from_file( key, filename );
}

void NCPA::StratifiedAtmosphere3D::get_minimum_altitude_limits( double &minlimit,
			double &maxlimit ) const {
	minlimit = profile_->get_minimum_altitude();
	maxlimit = minlimit;
}

void NCPA::StratifiedAtmosphere3D::get_maximum_altitude_limits( double &minlimit,
			double &maxlimit ) const {
	minlimit = profile_->get_maximum_altitude();
	maxlimit = minlimit;
}
