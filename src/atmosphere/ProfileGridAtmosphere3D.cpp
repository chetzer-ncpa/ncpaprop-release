#include "ProfileGridAtmosphere3D.h"
#include "AtmosphericProperty3D.h"
#include "AtmosphericProperty1D.h"
#include "Atmosphere1D.h"
#include "AtmosphericModel.h"
#include "LANLInterpolation.h"
#include "units.h"
#include "util.h"

#include <set>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include <cfloat>

#ifndef PI
#define PI std::acos( -1.0 )
#endif

// #ifndef GAMMA_FOR_C
// #define GAMMA_FOR_C 1.4
// #endif

// #ifndef R_FOR_C
// #define R_FOR_C 287.0
// #endif



NCPA::ProfileGridAtmosphere3D::ProfileGridAtmosphere3D( const std::string &summary_file,
	size_t nz, double *zvec, size_t precision, const std::string &header_file ) {
	read_atmosphere_from_file( summary_file, nz, zvec, precision, header_file );
}

NCPA::ProfileGridAtmosphere3D::ProfileGridAtmosphere3D( const std::string &summary_file,
	size_t nz, double *zvec ) {
	read_atmosphere_from_file( summary_file, nz, zvec, 0, "" );
}

NCPA::ProfileGridAtmosphere3D::ProfileGridAtmosphere3D( const std::string &summary_file,
	size_t nz, double *zvec, size_t precision ) {
	read_atmosphere_from_file( summary_file, nz, zvec, precision, "" );
}

NCPA::ProfileGridAtmosphere3D::ProfileGridAtmosphere3D( const std::string &summary_file,
	size_t nz, double *zvec, const std::string &header_file ) {
	read_atmosphere_from_file( summary_file, nz, zvec, 0, header_file );
}



NCPA::ProfileGridAtmosphere3D::~ProfileGridAtmosphere3D() {
	std::map< std::string, NCPA::VectorAtmosphericProperty3D * >::iterator it;
	for (it = vector_contents_.begin(); it != vector_contents_.end(); ++it) {
		delete it->second;
	}
	vector_contents_.clear();
	std::map< std::string, NCPA::ScalarAtmosphericProperty3D *>::iterator sit;
	for (sit = scalar_contents_.begin(); sit != scalar_contents_.end(); ++sit) {
		delete sit->second;
	}
	scalar_contents_.clear();
}


void NCPA::ProfileGridAtmosphere3D::read_atmosphere_from_file( const std::string &summary_file,
	size_t nz, double *zvec, size_t precision, const std::string &header_file ) {

	// read in the lines of the summary file
	std::ifstream infile( summary_file );
	std::string line;
	std::set<int> xvals, yvals;
	double scale = std::pow( 10.0, precision );
	std::vector< std::string > filenames, keys;
	std::vector<double> xvec, yvec;
	std::vector<std::string>::const_iterator cstr;
	std::map< int, int > xcounts;
	size_t i, j, nx, ny;
	bool zvec_from_profiles = (nz <= 0);

	NCPA::safe_getline( infile, line );
	while (line.size() > 0) {
		line = NCPA::deblank(line);
		if (line[0] != '#') {
			std::vector<std::string> fields = NCPA::split( line );
			if (fields.size() == 3) {
				//lines.push_back( line );
				xvec.push_back( std::stod( fields[ 0 ] ) );
				yvec.push_back( std::stod( fields[ 1 ] ) );
				filenames.push_back( fields[ 2 ] );
				int xval = (int)std::round( xvec.back() * scale );
				int yval = (int)std::round( yvec.back() * scale );
				// std::cout << "Inserting point [ " << xval << ", " << yval << " ]" << std::endl;
				xvals.insert( xval );
				yvals.insert( yval );
				xcounts[xval]++;
			}
		}
		NCPA::safe_getline( infile, line );
	}
	nx = xvals.size();
	ny = yvals.size();

	// step through xcounts and verify that all have the same number of y values associated
	int n_per_x = xcounts.cbegin()->second;
	std::map<int,int>::const_iterator cit = xcounts.cbegin();
	cit++;
	for ( ; cit != xcounts.cend(); ++cit) {
		if (cit->second != n_per_x) {
			throw std::runtime_error( "Inconsistent grid dimensions" );
		}
	}

	// Create x and y value arrays for storage and indexing
	double *x = new double[ nx ];
	std::memset( x, 0, nx * sizeof(double) );
	double *y = new double[ ny ];
	std::memset( y, 0, ny * sizeof(double) );
	i = 0;
	for (std::set<int>::const_iterator xit = xvals.cbegin(); xit != xvals.cend(); ++xit ) {
		// std::cout << "  x[ " << i << " ] = " << ((double)(*xit)) / scale << std::endl;
		x[ i++ ] = ((double)(*xit)) / scale;
	}
	i = 0;
	for (std::set<int>::const_iterator yit = yvals.cbegin(); yit != yvals.cend(); ++yit ) {
		// std::cout << "  y[ " << i << " ] = " << ((double)(*yit)) / scale << std::endl;
		y[ i++ ] = ((double)(*yit)) / scale;
	}

	// Create matrix of 1-D profiles
	NCPA::Atmosphere1D ***profilemat = NCPA::allocate_matrix<NCPA::Atmosphere1D *>( nx, ny );
	// for (cstr = lines.cbegin(); cstr != lines.cend(); ++cstr ) {
		// std::vector<std::string> fields = NCPA::split( line );
		// i = find_closest_index( x, nx, std::stod( fields[ 0 ] ) );
		// j = find_closest_index( y, ny, std::stod( fields[ 1 ] ) );
		// profilemat[ i ][ j ] = new NCPA::Atmosphere1D( fields[2], header_file );
	for (size_t lineind = 0; lineind < xvec.size(); lineind++) {

		i = NCPA::find_closest_index( x, nx, xvec[ lineind ] );
		j = NCPA::find_closest_index( y, ny, yvec[ lineind ] );
		profilemat[ i ][ j ] = new NCPA::Atmosphere1D( filenames[ lineind ], header_file );

		// make altitude units km internally for simplicity
		profilemat[ i ][ j ]->convert_altitude_units( NCPA::UNITS_DISTANCE_KILOMETERS );

		// std::cout << "Inserted atmosphere at [ " << i << ", " << j << " ]" << std::endl;
	}

	// if no z vector was specified, use the one at the origin
	if (zvec_from_profiles) {
		size_t xo = NCPA::find_closest_index( x, nx, 0.0 );
		size_t yo = NCPA::find_closest_index( y, ny, 0.0 );
		nz = profilemat[ xo ][ yo ]->nz();
		zvec = new double[ nz ];
		NCPA::units_t zunits;
		profilemat[ xo ][ yo ]->get_altitude_vector( zvec, &zunits );
		NCPA::Units::convert( zvec, nz, zunits, NCPA::UNITS_DISTANCE_KILOMETERS, zvec );
	}

	// create 2-D matrices of 1-D properties = 3-D properties
	// vectors first
	keys = profilemat[ 0 ][ 0 ]->get_vector_keys();
	NCPA::AtmosphericProperty1D ***propmat =
			NCPA::allocate_matrix<NCPA::AtmosphericProperty1D *>( nx, ny );
	for (cstr = keys.cbegin(); cstr != keys.cend(); ++cstr) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				propmat[ i ][ j ] = profilemat[ i ][ j ]->get_vector_property_object( *cstr );
			}
		}
		NCPA::VectorAtmosphericProperty3D *newprop = new NCPA::VectorAtmosphericProperty3D(
			*cstr, nx, x, ny, y, NCPA::UNITS_DISTANCE_KILOMETERS, propmat,
			nz, zvec, NCPA::UNITS_DISTANCE_KILOMETERS );
		add_property( *cstr, newprop );
	}
	NCPA::free_matrix( propmat, nx, ny );
	keys.clear();

	// now scalars
	keys = profilemat[ 0 ][ 0 ]->get_scalar_keys();
	NCPA::ScalarWithUnits ***spropmat =
			NCPA::allocate_matrix<NCPA::ScalarWithUnits *>( nx, ny );
	for (cstr = keys.cbegin(); cstr != keys.cend(); ++cstr) {
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				spropmat[ i ][ j ] = profilemat[ i ][ j ]->get_scalar_property_object( *cstr );
			}
		}
		NCPA::ScalarAtmosphericProperty3D *newsprop = new NCPA::ScalarAtmosphericProperty3D(
			*cstr, nx, x, ny, y, NCPA::UNITS_DISTANCE_KILOMETERS, spropmat );
		add_property( *cstr, newsprop );
	}
	NCPA::free_matrix( spropmat, nx, ny );
	keys.clear();

	NCPA::free_matrix_and_contents( profilemat, nx, ny );
	if (zvec_from_profiles) {
		delete [] zvec;
	}
}

void NCPA::ProfileGridAtmosphere3D::add_property( const std::string &key,
	NCPA::VectorAtmosphericProperty3D *prop ) {

	// make sure we don't already have it
	if (contains_property( key )) {
		throw new std::runtime_error( "Property " + key + " already exists!" );
	}
	std::string newkey( key );
	vector_contents_[ newkey ] = prop;
}

void NCPA::ProfileGridAtmosphere3D::add_property( const std::string &key,
	NCPA::ScalarAtmosphericProperty3D *prop ) {

	// make sure we don't already have it
	if (contains_property( key )) {
		throw new std::runtime_error( "Property " + key + " already exists!" );
	}
	std::string newkey( key );
	scalar_contents_[ newkey ] = prop;
}

void NCPA::ProfileGridAtmosphere3D::convert_range_units( NCPA::units_t new_units ) {

	std::map<std::string, NCPA::VectorAtmosphericProperty3D *>::iterator vit;
	std::map<std::string, NCPA::ScalarAtmosphericProperty3D *>::iterator sit;
	for (vit = vector_contents_.begin(); vit != vector_contents_.end(); ++vit) {
		vit->second->convert_range_units( new_units );
	}
	for (sit = scalar_contents_.begin(); sit != scalar_contents_.end(); ++sit) {
		sit->second->convert_range_units( new_units );
	}
}

void NCPA::ProfileGridAtmosphere3D::convert_altitude_units( NCPA::units_t new_units ) {

	std::map<std::string, NCPA::VectorAtmosphericProperty3D *>::iterator vit;
	for (vit = vector_contents_.begin(); vit != vector_contents_.end(); ++vit) {
		vit->second->convert_altitude_units( new_units );
	}
}

void NCPA::ProfileGridAtmosphere3D::convert_property_units( const std::string &key,
	NCPA::units_t new_units ) {

	std::map<std::string, NCPA::VectorAtmosphericProperty3D *>::iterator vit;
	std::map<std::string, NCPA::ScalarAtmosphericProperty3D *>::iterator sit;
	vit = vector_contents_.find( key );
	if (vit != vector_contents_.end()) {
		vit->second->convert_property_units( new_units );
	} else {
		sit = scalar_contents_.find( key );
		if (sit != scalar_contents_.end()) {
			sit->second->convert_property_units( new_units );
		} else {
			throw std::runtime_error( "Requested key " + key + " not found" );
		}
	}
}

bool NCPA::ProfileGridAtmosphere3D::contains_vector( const std::string &key ) const {
	std::map<std::string, NCPA::VectorAtmosphericProperty3D *>::const_iterator vit
		= vector_contents_.find( key );
	return (vit != vector_contents_.end());
}

bool NCPA::ProfileGridAtmosphere3D::contains_scalar( const std::string &key ) const {
	std::map<std::string, NCPA::ScalarAtmosphericProperty3D *>::const_iterator sit
		= scalar_contents_.find( key );
	return (sit != scalar_contents_.end());
}

bool NCPA::ProfileGridAtmosphere3D::contains_property( const std::string &key ) const {
	return contains_vector( key ) || contains_scalar( key );
}

double NCPA::ProfileGridAtmosphere3D::get( double x, double y, const std::string &key ) const {
	std::map< std::string, NCPA::ScalarAtmosphericProperty3D *>::const_iterator
			it = scalar_contents_.find( key );
	if (it == scalar_contents_.cend()) {
		throw new std::runtime_error( "Key " + key + " not found." );
	}

	return it->second->get( x, y );
}

double NCPA::ProfileGridAtmosphere3D::get( double x, double y, double z,
			const std::string &key ) const {
	std::map< std::string, NCPA::VectorAtmosphericProperty3D *>::const_iterator
			it = vector_contents_.find( key );
	if (it == vector_contents_.cend()) {
		throw new std::runtime_error( "Key " + key + " not found." );
	}

	return it->second->get( x, y, z );
}

NCPA::units_t NCPA::ProfileGridAtmosphere3D::get_range_units() const {
	return vector_contents_.cbegin()->second->get_range_units();
}

NCPA::units_t NCPA::ProfileGridAtmosphere3D::get_altitude_units() const {
	return vector_contents_.cbegin()->second->get_altitude_units();
}

NCPA::units_t NCPA::ProfileGridAtmosphere3D::get_property_units(
	const std::string &key ) const {

	if (contains_vector( key )) {
		return vector_contents_.find( key )->second->get_property_units();
	} else if (contains_scalar( key )) {
		return scalar_contents_.find( key )->second->get_property_units();
	} else {
		return UNITS_NONE;
	}
}

double NCPA::ProfileGridAtmosphere3D::get_maximum_valid_range( double azimuth ) const {
	// std::cout << "Azimuth: " << azimuth << std::endl;
	double az_math = NCPA::deg2rad( 90.0 - azimuth );
	while (az_math < 0) {
		az_math += 2.0 * PI;
	}
	while (az_math >= (2.0 * PI)) {
		az_math -= 2.0 * PI;
	}

	// get horizontal limits
	double xmax = vector_contents_.cbegin()->second->x_max();
	double xmin = vector_contents_.cbegin()->second->x_min();
	double ymax = vector_contents_.cbegin()->second->y_max();
	double ymin = vector_contents_.cbegin()->second->y_min();
	// std::cout << "Limits: " << xmin << ", " << xmax << ", "
	// 		  << ymin << ", " << ymax << std::endl;

	// quadrant points
	double quadrant[ 4 ];
	quadrant[ 0 ] = NCPA::normalizeAzimuthRadians( std::atan( ymax / xmax ) );
	quadrant[ 1 ] = NCPA::normalizeAzimuthRadians( PI - std::atan( ymax / -xmin ) );
	quadrant[ 2 ] = NCPA::normalizeAzimuthRadians( PI + std::atan( ymin / xmin ) );
	quadrant[ 3 ] = NCPA::normalizeAzimuthRadians( -std::atan( -ymin / xmax ) );

	if (az_math <= quadrant[0] || az_math > quadrant[ 3 ]) {
		return xmax / std::cos( az_math );
	} else if (az_math > quadrant[ 0 ] && az_math <= quadrant[ 1 ]) {
		return ymax / std::sin( az_math );
	} else if (az_math > quadrant[ 1 ] && az_math <= quadrant[ 2 ]) {
		return xmin / std::cos( az_math );
	} else {
		return ymin / std::sin( az_math );
	}
}

double NCPA::ProfileGridAtmosphere3D::get_derivative( double x, double y,
	const std::string &key, size_t order, deriv_t *directions ) const {

	std::map< std::string, NCPA::ScalarAtmosphericProperty3D * >::const_iterator cit;
	cit = scalar_contents_.find( key );
	if (cit != scalar_contents_.cend()) {
		return cit->second->get_derivative( x, y, order, directions );
	} else {
		throw std::runtime_error( "Scalar key " + key + " not found" );
	}
}

double NCPA::ProfileGridAtmosphere3D::get_derivative( double x, double y, double z,
			const std::string &key, size_t order, deriv_t *directions ) const {

	std::map< std::string, NCPA::VectorAtmosphericProperty3D * >::const_iterator cit;
	cit = vector_contents_.find( key );
	if (cit != vector_contents_.cend()) {
		return cit->second->get_derivative( x, y, z, order, directions );
	} else {
		throw std::runtime_error( "Vector key " + key + " not found" );
	}
}

// add 3-D array of values, assuming they're at the same grid points as
void NCPA::ProfileGridAtmosphere3D::add_property( const std::string &key, double ***prop,
	size_t nxi, size_t nyi, size_t nzi, NCPA::units_t prop_units ) {

	if (contains_property(key)) {
		throw std::runtime_error( "Key " + key + " already exists" );
	}

	NCPA::VectorAtmosphericProperty3D *basis = vector_contents_.begin()->second;
	double *x, *y, *z;
	size_t nx, ny, nz;
	basis->x_vector( nx, x );
	basis->y_vector( ny, y );
	basis->z_vector( nz, z );

	// check for dimension agreement
	if (nx != nxi || ny != nyi || nz != nzi) {
		throw std::runtime_error( "Dimensional conflict in parameter " + key );
	}

	vector_contents_[ key ] = new NCPA::VectorAtmosphericProperty3D(
		key, nx, x, ny, y, nz, z, prop, basis->get_range_units(),
		basis->get_altitude_units(), prop_units );
	delete [] x;
	delete [] y;
	delete [] z;
}


void NCPA::ProfileGridAtmosphere3D::add_property( const std::string &key, double **prop,
	size_t nxi, size_t nyi, NCPA::units_t units ) {

	if (contains_property(key)) {
		throw std::runtime_error( "Key " + key + " already exists" );
	}

	NCPA::ScalarAtmosphericProperty3D *basis = scalar_contents_.begin()->second;
	double *x, *y;
	size_t nx, ny;
	basis->x_vector( nx, x );
	basis->y_vector( ny, y );
	// check for dimension agreement
	if (nx != nxi || ny != nyi) {
		throw std::runtime_error( "Dimensional conflict in parameter " + key );
	}

	scalar_contents_[ key ] = new NCPA::ScalarAtmosphericProperty3D(
		key, nx, x, ny, y, prop, basis->get_range_units(), units );
	delete [] x;
	delete [] y;
}


void NCPA::ProfileGridAtmosphere3D::calculate_sound_speed_from_temperature(
	const std::string &new_key,
	const std::string &temperature_key, NCPA::units_t wind_units ) {

	if (!contains_vector(temperature_key)) {
		throw std::runtime_error( "No property for temperature key "
			+ temperature_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator t_it;
	// index variables
	size_t nx, ny, nz, i, j, k;
	double ***t, ***c;

	// compute sound speed and convert to desired units
	t_it = vector_contents_.find(temperature_key);
	t_it->second->as_matrix( t, nx, ny, nz );
	NCPA::units_t temperature_units = t_it->second->get_property_units();
	c = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				c[ i ][ j ][ k ] = NCPA::Units::convert(
					NCPA::AtmosphericModel::soundspeed_from_temperature(
						NCPA::Units::convert(
							t[ i ][ j ][ k ], temperature_units,
							NCPA::UNITS_TEMPERATURE_KELVIN
						) ),
					NCPA::UNITS_SPEED_METERS_PER_SECOND, wind_units );
			}
		}
	}

	add_property( new_key, c, nx, ny, nz, wind_units );
	NCPA::free_matrix3d<double>( c, nx, ny, nz );
	t_it->second->free_matrix( t );
}


void NCPA::ProfileGridAtmosphere3D::calculate_sound_speed_from_pressure_and_density(
	const std::string &new_key, const std::string &pressure_key,
	const std::string &density_key, units_t wind_units ) {

	if (!contains_vector(pressure_key)) {
		throw std::runtime_error( "No property for pressure key "
			+ pressure_key + " found." );
	}
	if (!contains_vector(density_key)) {
		throw std::runtime_error( "No property for density key "
			+ density_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}

	// input variables
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator p_it, d_it;
	size_t nx, ny, nz, i, j, k;
	double ***p, ***d, ***c;
	NCPA::units_t p_units, d_units;
	p_it = vector_contents_.find( pressure_key );
	p_it->second->as_matrix( p, nx, ny, nz );
	p_units = p_it->second->get_property_units();
	d_it = vector_contents_.find( density_key );
	d_it->second->as_matrix( d, nx, ny, nz );
	d_units = d_it->second->get_property_units();

	c = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				c[ i ][ j ][ k ] = NCPA::Units::convert(
					NCPA::AtmosphericModel::soundspeed_from_pressure_density(
						NCPA::Units::convert(
							p[ i ][ j ][ k ], p_units,
							NCPA::UNITS_PRESSURE_PASCALS ),
						NCPA::Units::convert(
							d[ i ][ j ][ k ], d_units,
							NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER)
						),
					NCPA::UNITS_SPEED_METERS_PER_SECOND, wind_units );
			}
		}
	}

	add_property( new_key, c, nx, ny, nz, wind_units );
	NCPA::free_matrix3d<double>( c, nx, ny, nz );
	p_it->second->free_matrix( p );
	d_it->second->free_matrix( d );
}


void NCPA::ProfileGridAtmosphere3D::calculate_wind_speed(
	const std::string &new_key,
	const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key ) {

	if (!contains_vector(we_wind_speed_key)) {
		throw std::runtime_error( "No property for zonal wind key "
			+ we_wind_speed_key + " found." );
	}
	if (!contains_vector(sn_wind_speed_key)) {
		throw std::runtime_error( "No property for meridional wind key "
			+ sn_wind_speed_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}

	// input variables
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator u_it, v_it;
	size_t nx, ny, nz, i, j, k;
	double ***u, ***v, ***c;
	NCPA::units_t u_units, v_units;
	u_it = vector_contents_.find( we_wind_speed_key );
	u_it->second->as_matrix( u, nx, ny, nz );
	u_units = u_it->second->get_property_units();
	v_it = vector_contents_.find( sn_wind_speed_key );
	v_it->second->as_matrix( v, nx, ny, nz );
	v_units = v_it->second->get_property_units();
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed components" );
	}

	c = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				c[ i ][ j ][ k ] = std::sqrt(
					u[ i ][ j ][ k ] * u[ i ][ j ][ k ]
					+ v[ i ][ j ][ k ] * v[ i ][ j ][ k ] );
			}
		}
	}

	add_property( new_key, c, nx, ny, nz, u_units );
	NCPA::free_matrix3d<double>( c, nx, ny, nz );
	u_it->second->free_matrix( u );
	v_it->second->free_matrix( v );
}


void NCPA::ProfileGridAtmosphere3D::calculate_wind_direction(
	const std::string &new_key,
	const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
	units_t direction_units ) {

	if (!contains_vector(we_wind_speed_key)) {
		throw std::runtime_error( "No property for zonal wind key "
			+ we_wind_speed_key + " found." );
	}
	if (!contains_vector(sn_wind_speed_key)) {
		throw std::runtime_error( "No property for meridional wind key "
			+ sn_wind_speed_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}

	// input variables
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator u_it, v_it;
	size_t nx, ny, nz, i, j, k;
	double ***u, ***v, ***c;
	NCPA::units_t u_units, v_units;
	u_it = vector_contents_.find( we_wind_speed_key );
	u_it->second->as_matrix( u, nx, ny, nz );
	u_units = u_it->second->get_property_units();
	v_it = vector_contents_.find( sn_wind_speed_key );
	v_it->second->as_matrix( v, nx, ny, nz );
	v_units = v_it->second->get_property_units();
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed components" );
	}

	c = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				c[ i ][ j ][ k ] =
					NCPA::Units::convert(
						NCPA::Units::convert(
							PI/2.0 - std::atan2( v[ i ][ j ][ k ], u[ i ][ j ][ k ] ),
							NCPA::UNITS_ANGLE_RADIANS, NCPA::UNITS_ANGLE_DEGREES
						),
						NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, direction_units );
			}
		}
	}
	add_property( new_key, c, nx, ny, nz, direction_units );
	NCPA::free_matrix3d<double>( c, nx, ny, nz );
	u_it->second->free_matrix( u );
	v_it->second->free_matrix( v );
}


void NCPA::ProfileGridAtmosphere3D::calculate_wind_component(
	const std::string &new_key,
	const std::string &wind_speed_key, const std::string &wind_direction_key,
	double azimuth ) {

	if (!contains_vector(wind_speed_key)) {
		throw std::runtime_error( "No property for wind speed key "
			+ wind_speed_key + " found." );
	}
	if (!contains_vector(wind_direction_key)) {
		throw std::runtime_error( "No property for wind direction key "
			+ wind_direction_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}

	// input variables
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator s_it, d_it;
	size_t nx, ny, nz, i, j, k;
	double ***s, ***d, ***c;
	NCPA::units_t s_units, d_units;
	s_it = vector_contents_.find( wind_speed_key );
	s_it->second->as_matrix( s, nx, ny, nz );
	s_units = s_it->second->get_property_units();
	d_it = vector_contents_.find( wind_direction_key );
	d_it->second->as_matrix( d, nx, ny, nz );
	d_units = d_it->second->get_property_units();
	double az_rad = NCPA::Units::convert( azimuth,
		NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS );

	c = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				c[ i ][ j ][ k ] =
					s[ i ][ j ][ k ] * std::cos(
						NCPA::Units::convert( d[ i ][ j ][ k ],
							NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS
						) - az_rad
					);
			}
		}
	}
	add_property( new_key, c, nx, ny, nz, s_units );
	NCPA::free_matrix3d<double>( c, nx, ny, nz );
	s_it->second->free_matrix( s );
	d_it->second->free_matrix( d );
}


void NCPA::ProfileGridAtmosphere3D::calculate_effective_sound_speed(
	const std::string &new_key,
	const std::string &sound_speed_key, const std::string &wind_component_key ) {

	if (!contains_vector(sound_speed_key)) {
		throw std::runtime_error( "No property for sound speed key "
			+ sound_speed_key + " found." );
	}
	if (!contains_vector(wind_component_key)) {
		throw std::runtime_error( "No property for wind component key "
			+ wind_component_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}

	// input variables
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator s_it, w_it;
	size_t nx, ny, nz, i, j, k;
	double ***s, ***w, ***c;
	NCPA::units_t s_units, w_units;
	s_it = vector_contents_.find( sound_speed_key );
	s_it->second->as_matrix( s, nx, ny, nz );
	s_units = s_it->second->get_property_units();
	w_it = vector_contents_.find( wind_component_key );
	w_it->second->as_matrix( w, nx, ny, nz );
	w_units = w_it->second->get_property_units();

	c = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				// use convert in case of units mismatch
				c[ i ][ j ][ k ] =
					s[ i ][ j ][ k ] + NCPA::Units::convert(
						w[ i ][ j ][ k ], w_units, s_units );
			}
		}
	}
	add_property( new_key, c, nx, ny, nz, s_units );
	NCPA::free_matrix3d<double>( c, nx, ny, nz );
	s_it->second->free_matrix( s );
	w_it->second->free_matrix( w );
}


void NCPA::ProfileGridAtmosphere3D::calculate_attenuation(
	const std::string &new_key,
	const std::string &temperature_key, const std::string &pressure_key,
	const std::string &density_key, double freq, double tweak_factor ) {

	if (!contains_vector(temperature_key)) {
		throw std::runtime_error( "No property for temperature key "
			+ temperature_key + " found." );
	}
	if (!contains_vector(pressure_key)) {
		throw std::runtime_error( "No property for pressure key "
			+ pressure_key + " found." );
	}
	if (!contains_vector(density_key)) {
		throw std::runtime_error( "No property for density key "
			+ density_key + " found." );
	}
	if (contains_property(new_key)) {
		throw std::runtime_error( "Property with key " + new_key + " already exists" );
	}

	// input variables
	std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::iterator t_it, p_it, d_it;
	size_t nx, ny, nz, i, j, k;
	double ***t, ***p, ***d, ***a, *z;
	NCPA::units_t t_units, p_units, d_units, z_units;
	t_it = vector_contents_.find( temperature_key );
	t_it->second->as_matrix( t, nx, ny, nz );
	t_units = t_it->second->get_property_units();
	p_it = vector_contents_.find( pressure_key );
	p_it->second->as_matrix( p, nx, ny, nz );
	p_units = p_it->second->get_property_units();
	d_it = vector_contents_.find( density_key );
	d_it->second->as_matrix( d, nx, ny, nz );
	d_units = d_it->second->get_property_units();
	t_it->second->z_vector( nz, z );
	z_units = t_it->second->get_altitude_units();

	a = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				// use convert in case of units mismatch
				a[ i ][ j ][ k ] =
					NCPA::AtmosphericModel::attenuation_sutherland_bass(
						NCPA::Units::convert( z[ k ],
							z_units, NCPA::UNITS_DISTANCE_KILOMETERS ),
						NCPA::Units::convert( t[ i ][ j ][ k ],
							t_units, NCPA::UNITS_TEMPERATURE_KELVIN ),
						NCPA::Units::convert( p[ i ][ j ][ k ],
							p_units, NCPA::UNITS_PRESSURE_PASCALS ),
						NCPA::Units::convert( d[ i ][ j ][ k ],
							d_units, NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER ),
						freq ) * tweak_factor;
			}
		}
	}
	add_property( new_key, a, nx, ny, nz, NCPA::UNITS_NONE );
	NCPA::free_matrix3d<double>( a, nx, ny, nz );
	t_it->second->free_matrix( t );
	p_it->second->free_matrix( p );
	d_it->second->free_matrix( d );
	delete [] z;
}


void NCPA::ProfileGridAtmosphere3D::remove_property( const std::string &key ) {
	remove_vector_property( key );
	remove_scalar_property( key );
}

void NCPA::ProfileGridAtmosphere3D::remove_vector_property( const std::string &key ) {
	if (contains_vector( key ) ) {
		std::map<std::string,NCPA::VectorAtmosphericProperty3D *>::iterator it
			= vector_contents_.find( key );
		delete it->second;
		vector_contents_.erase( it );
	}
}

void NCPA::ProfileGridAtmosphere3D::remove_scalar_property( const std::string &key ) {
	if (contains_scalar( key ) ) {
		std::map<std::string,NCPA::ScalarAtmosphericProperty3D *>::iterator it
			= scalar_contents_.find( key );
		delete it->second;
		scalar_contents_.erase( it );
	}
}

double NCPA::ProfileGridAtmosphere3D::get_minimum_altitude( double x, double y ) const {
	double minz = DBL_MAX;
	for (std::map<std::string,NCPA::VectorAtmosphericProperty3D*>::const_iterator
		cit = vector_contents_.begin(); cit != vector_contents_.cend(); ++cit) {
		minz = NCPA::min( minz, cit->second->z_min() );
	}
	return minz;
}


void NCPA::ProfileGridAtmosphere3D::get_property_template( const std::string &basis,
				size_t &nx, double *x, NCPA::units_t &x_units,
				size_t &ny, double *y, NCPA::units_t &y_units,
				double **&prop ) const {

	NCPA::AtmosphericProperty3D *propptr = get_property( basis );

	propptr->x_vector( nx, x );
	x_units = propptr->get_range_units();
	propptr->y_vector( ny, y );
	y_units = x_units;
	prop = NCPA::allocate_matrix<double>( nx, ny );
}

void NCPA::ProfileGridAtmosphere3D::get_property_template( const std::string &basis,
			size_t &nx, double *x, NCPA::units_t &x_units,
			size_t &ny, double *y, NCPA::units_t &y_units,
			size_t &nz, double *z, NCPA::units_t &z_units,
			double ***&prop ) const {
	NCPA::AtmosphericProperty3D *propptr = get_property( basis );

	propptr->x_vector( nx, x );
	x_units = propptr->get_range_units();
	propptr->y_vector( ny, y );
	y_units = x_units;
	propptr->z_vector( nz, z );
	z_units = propptr->get_altitude_units();
	prop = NCPA::matrix3d<double>( nx, ny, nz );

}

void NCPA::ProfileGridAtmosphere3D::free_property_template(
				size_t nx, double *x,
				size_t ny, double *y,
				double **prop ) const {
	NCPA::free_matrix<double>( prop, nx, ny );
	delete [] x;
	delete [] y;
}

void NCPA::ProfileGridAtmosphere3D::free_property_template(
				size_t nx, double *x,
				size_t ny, double *y,
				size_t nz, double *z,
				double ***prop ) const {
	NCPA::free_matrix3d<double>( prop, nx, ny, nz );
	delete [] x;
	delete [] y;
	delete [] z;
}

NCPA::AtmosphericProperty3D *NCPA::ProfileGridAtmosphere3D::get_property(
	const std::string &key ) const {
	if (contains_vector(key)) {
		return (NCPA::AtmosphericProperty3D *)(vector_contents_.find(key)->second);
	} else if (contains_scalar(key)) {
		return (NCPA::AtmosphericProperty3D *)(scalar_contents_.find(key)->second);
	}
	throw std::runtime_error( "No property with key " + key + " found." );
}

void NCPA::ProfileGridAtmosphere3D::copy_property( const std::string &old_key,
				const std::string &new_key ) {
	if (contains_property(new_key)) {
		throw std::runtime_error( "Requested new key " + new_key + " already exists" );
	}
	if (contains_vector(old_key)) {
		vector_contents_[new_key] = new NCPA::VectorAtmosphericProperty3D(
			*(vector_contents_.find(old_key)->second) );
	} else if (contains_scalar(old_key)) {
		scalar_contents_[new_key] = new NCPA::ScalarAtmosphericProperty3D(
			*(scalar_contents_.find(old_key)->second) );
	}
	throw std::runtime_error( "Key " + old_key + " not found" );
}

void NCPA::ProfileGridAtmosphere3D::read_attenuation_from_file( const std::string &key,
			const std::string &filename ) {

	std::ifstream in( filename );
	std::string line;
	std::ostringstream oss;
	std::getline( in, line );
	std::vector< std::string > atmlines, fields;
	size_t i, j, k;

	if (contains_vector( key )) {
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	}

	while ( in.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = NCPA::deblank( line );
		if (line[ 0 ] != '#') {
			atmlines.push_back( line );
		}

		getline( in, line );
	}
	in.close();

	size_t nlines = atmlines.size();
	double *z_a = new double[ nlines ];
	std::memset( z_a, 0, nlines * sizeof( double ) );
	double *attn = new double[ nlines ];
	std::memset( attn, 0, nlines * sizeof( double ) );
	double this_z, this_a;

	for (i = 0; i < nlines; i++) {
		fields = NCPA::split( NCPA::deblank( atmlines[ i ] ), " ," );
		if (fields.size() != 2) {
			oss << "ProfileGridAtmosphere3D - Error parsing attenuation line:" << std::endl << line << std::endl
				<< "Must be formatted as:" << std::endl
				<< "altitude  attenuation" << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		try {
			this_z = std::stof( fields[ 0 ] );
			this_a = std::stof( fields[ 1 ] );
		} catch ( std::invalid_argument &e ) {
			oss << "Atmosphere1D - Error parsing attenuation line:" << std::endl << line << std::endl
				<< "Both fields must be numerical" << std::endl;
			throw std::invalid_argument( oss.str() );
		}
		z_a[ i ] = this_z;
		attn[ i ] = this_a;
	}

	// get new property template
	double ***a_mat, *x_vec, *y_vec, *z_vec;
	size_t nx_vec, ny_vec, nz_vec;
	NCPA::units_t xunits, yunits, z_units, alt_units;
	alt_units = get_altitude_units();
	convert_altitude_units( NCPA::Units::fromString("km" ) );
	get_property_template( "T",
		nx_vec, x_vec, xunits, ny_vec, y_vec, yunits, nz_vec, z_vec, z_units, a_mat );

	// extend if necessary
	bool short_on_bottom = (z_a[0] > z_vec[0]);
	bool short_on_top    = (z_a[nlines-1] < z_vec[ nz_vec-1 ]);
	double *temp_a, *temp_z;
	if (short_on_bottom) {
		temp_z = new double[ nlines+1 ];
		temp_a = new double[ nlines+1 ];
		std::memcpy( temp_z+1, z_a, nlines );
		std::memcpy( temp_a+1, attn, nlines );
		temp_z[ 0 ] = z_vec[ 0 ];
		temp_a[ 0 ] = temp_a[ 1 ];
		delete [] z_a;
		delete [] attn;
		z_a = temp_z;
		attn = temp_a;
		nlines++;
	}
	if (short_on_top) {
		temp_z = new double[ nlines+1 ];
		temp_a = new double[ nlines+1 ];
		std::memcpy( temp_z, z_a, nlines );
		std::memcpy( temp_a, attn, nlines );
		temp_z[ nlines-1 ] = z_vec[ nz_vec-1 ];
		temp_a[ nlines-1 ] = temp_a[ nlines-2 ];
		delete [] z_a;
		delete [] attn;
		z_a = temp_z;
		attn = temp_a;
		nlines++;
	}

	// Now interpolate onto current z grid
	LANL::natural_cubic_spline_1D aspl;
	LANL::prep( aspl, nlines );
	std::memcpy( aspl.x_vals, z_a, nlines*sizeof(double) );
	std::memcpy( aspl.f_vals, attn, nlines*sizeof(double) );
	LANL::set( aspl );
	for (k = 0; k < nz_vec; k++) {
		double attn_z = LANL::eval_f( z_vec[ k ], aspl );
		for (i = 0; i < nx_vec; i++) {
			for (j = 0; j < ny_vec; j++) {
				a_mat[ i ][ j ][ k ] = attn_z;
			}
		}
	}
	LANL::clear( aspl );

	add_property( key, a_mat, nx_vec, ny_vec, nz_vec, NCPA::UNITS_NONE );
	free_property_template( nx_vec, x_vec, ny_vec, y_vec, nz_vec, z_vec, a_mat );
	convert_altitude_units( alt_units );
}

void NCPA::ProfileGridAtmosphere3D::get_minimum_altitude_limits( double &minvalue,
			double &maxvalue ) const {
	std::map<std::string,NCPA::VectorAtmosphericProperty3D *>::const_iterator it
			= vector_contents_.cbegin();
	minvalue = it->second->z_min();
	maxvalue = minvalue;
	++it;

	for ( ; it != vector_contents_.cend(); ++it) {
		minvalue = NCPA::min( minvalue, it->second->z_min() );
		maxvalue = NCPA::max( maxvalue, it->second->z_min() );
	}
}

void NCPA::ProfileGridAtmosphere3D::get_maximum_altitude_limits( double &minvalue,
			double &maxvalue ) const {
	std::map<std::string,NCPA::VectorAtmosphericProperty3D *>::const_iterator it
			= vector_contents_.cbegin();
	minvalue = it->second->z_max();
	maxvalue = minvalue;
	++it;

	for ( ; it != vector_contents_.cend(); ++it) {
		minvalue = NCPA::min( minvalue, it->second->z_max() );
		maxvalue = NCPA::max( maxvalue, it->second->z_max() );
	}
}
