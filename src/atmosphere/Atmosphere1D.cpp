//#include "Atmosphere.h"
#include "Atmosphere1D.h"
#include "AtmosphericProperty1D.h"
#include "units.h"
#include "util.h"
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

// #ifndef GAMMA_FOR_C
// #define GAMMA_FOR_C 1.4
// #endif

// #ifndef R_FOR_C
// #define R_FOR_C 287.0
// #endif

#ifndef PI
#define PI 3.14159
#endif

NCPA::Atmosphere1D::Atmosphere1D() {
	contents_.clear();
	scalar_contents_.clear();
	z_ = NULL;
}

NCPA::Atmosphere1D::Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units ) {
	contents_.clear();
	scalar_contents_.clear();
	z_ = new NCPA::VectorWithUnits( n_altitude_points, altitude_points, altitude_units );

	//nz_ = n_altitude_points;
	//z_ = new double[ nz_ ];
	//std::memcpy( z_, altitude_points, nz_ * sizeof(double) );
	//z_units_.push( altitude_units );
}

// NCPA::Atmosphere1D::Atmosphere1D( std::istream& in ) {
// 	read_from_stream( in );
// }

NCPA::Atmosphere1D::Atmosphere1D( const std::string &filename,
	const std::string &headerfilename ) {

	std::string hfile( headerfilename );
	if (headerfilename.size() == 0) {
		hfile = filename;
	}
	std::ifstream header_in( hfile );
	read_header_from_stream( header_in );
	header_in.close();

	std::ifstream in( filename );
	read_values_from_stream( in );
	in.close();

	headerlines.clear();
}

NCPA::Atmosphere1D::Atmosphere1D( const Atmosphere1D &source ) {
	z_ = new NCPA::VectorWithUnits( *(source.z_) );
	for ( std::map< std::string, NCPA::AtmosphericProperty1D * >::const_iterator it = source.contents_.cbegin();
		  it != source.contents_.cend(); ++it ) {
		NCPA::AtmosphericProperty1D *new_prop = new NCPA::AtmosphericProperty1D( *(it->second) );
		contents_[ it->first ] = new_prop;
	}
	for ( std::map< std::string, NCPA::ScalarWithUnits * >::const_iterator it = source.scalar_contents_.cbegin();
		  it != source.scalar_contents_.cend(); ++it ) {
		NCPA::ScalarWithUnits *new_scalar = new NCPA::ScalarWithUnits( *(it->second) );
		scalar_contents_[ it->first ] = new_scalar;
	}
}

void NCPA::Atmosphere1D::read_header_from_stream( std::istream& in ) {

	std::string line;

	std::getline( in, line );
	while ( in.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = NCPA::deblank( line );
		if (line[ 0 ] == '#') {
			// check second character
			if (line.size() > 1 && line[ 1 ] == '%') {
				headerlines.push_back( line.substr( 2 ) );
			} // otherwise it's a regular comment and can be ignored

		}

		std::getline( in, line );
	}

}

void NCPA::Atmosphere1D::read_values_from_stream( std::istream& in ) {
	if ( ! in.good() ) {
		throw std::runtime_error( "Atmosphere1D - Input stream not in good state" );
	}

	std::string line;
	std::vector< std::string > atmlines, scalarlines;
	std::ostringstream oss;     // for exceptions
	size_t i;                   // repeated index variable

	std::getline( in, line );
	while ( in.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = NCPA::deblank( line );
		if (line[ 0 ] == '#') {
			// check second character
			// if (line.size() > 1 && line[ 1 ] == '%') {
			// 	headerlines.push_back( line.substr( 2 ) );
			// } // otherwise it's a regular comment and can be ignored

		} else if (line.size() == 0) {
			// skip empty lines
		} else {
			atmlines.push_back( line );
		}

		std::getline( in, line );
	}
	//in.close();
	//cout << "Found " << headerlines.size() << " header lines" << endl;
	//cout << "Found " << atmlines.size() << " data lines" << endl;

	// parse them out
	size_t nfields = headerlines.size();
	if (nfields == 0) {
		throw std::runtime_error( "Atmosphere1D - No descriptive fields found." );
	}

	// hold contents
	std::vector< std::string > keys, fields;
	std::vector< unsigned int > column_numbers;
	std::vector< double > values;
	std::vector< NCPA::units_t > units;

	for (std::vector< std::string >::const_iterator it = headerlines.begin(); it != headerlines.end(); ++it) {
		fields = NCPA::split( NCPA::deblank( *it ), " ," );
		if ( fields.size() != 3 && fields.size() != 4 ) {
			oss << "Atmosphere1D - Error parsing descriptive line:" << std::endl << line << std::endl 
				<< "Must be formatted as:" << std::endl
				<< "column,key,units[,value]" << std::endl
				<< "Use column=0 and specify value for scalar quantities." << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		// process fields
		// field line needs to be parseable as an integer
		unsigned int col;
		try {
			col = (unsigned int)std::stoi( fields[ 0 ] );
		} catch ( std::invalid_argument &e ) {
			oss << "Atmosphere1D - Error parsing descriptive line:" << std::endl << line << std::endl 
				<< "First field not parseable as an integer" << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		double tempval = 0.0;
		//NCPA::units_t tempunits = NCPA::UNITS_NONE;
		NCPA::units_t tempunits = NCPA::Units::fromString( NCPA::deblank( fields[ 2 ] ) );
		if (tempunits == NCPA::UNITS_NONE) {
			throw std::invalid_argument( "Unrecognized units string: " + fields[ 2 ] );
		}
/*
		try {
			tempunits = NCPA::Units::fromString( NCPA::deblank( fields[ 2 ] ) );
		} catch (std::out_of_range& oor) {
			throw std::invalid_argument( oor.what() );
		}
*/
		if (fields.size() == 4) {
			tempval = std::stof( deblank(fields[ 3 ]) );
		}

		// add to header vectors
		column_numbers.push_back( col );
		keys.push_back( deblank( fields[ 1 ] ) );
		units.push_back( tempunits );
		values.push_back( tempval );  // this will be ignored for vector quantities

		//cout << "Column " << col << ": " << deblank(fields[1]) << " [ " << Units::toStr( tempunits ) << " ] = " << tempval << endl;
	}

	// check for uniqueness in keys
	std::vector< std::string > tempkeys( keys );
	std::sort( tempkeys.begin(), tempkeys.end() );
	std::vector< std::string >::iterator uit = std::unique( tempkeys.begin(), tempkeys.end() );
	if (uit != tempkeys.end()) {
		throw std::invalid_argument( "Atmosphere1D: key values are not unique" );
	}


	int nlines = atmlines.size();
	unsigned int ncols = 0;
	NCPA::units_t depunits = NCPA::UNITS_NONE;
	for (i = 0; i < column_numbers.size(); i++) {
		ncols = column_numbers[ i ] > ncols ? column_numbers[ i ] : ncols;
		if (column_numbers[ i ] == 1) {
			depunits = units[ i ];
		}
	}

	std::vector< double * > columns( ncols );
	for ( i = 0; i < ncols; i++ ) {
		double *col = new double[ nlines ];
		columns[ i ] = col;
	}


	// step through the data lines
	size_t row = 0;
	std::vector< std::string >::const_iterator it;
	for (it = atmlines.cbegin(); it != atmlines.cend(); ++it ) {
		fields.clear();
		fields = NCPA::split( NCPA::deblank( *it ), " \t," );
		if (fields.size() != ncols ) {
			oss << "Error parsing data line:" << std::endl << *it << std::endl
				<< ncols << " columns expected, " << fields.size() << " columns found.";
			for ( i = 0; i < ncols; i++ ) {
				delete [] columns[ i ];
			}
			throw std::invalid_argument( oss.str() );
		}
		for ( i = 0; i < ncols; i++ ) {
			try {
				double *thiscol = columns[ i ];
				thiscol[ row ] = std::stof( fields[ i ] );
			} catch (std::invalid_argument &e) {
				oss << "Error parsing data line:" << std::endl << *it << std::endl
					<< "Can't parse field " << fields[ i ] << " as a double";
				for ( i = 0; i < ncols; i++ ) {
					delete [] columns[ i ];
				}
				throw std::invalid_argument( oss.str() );
			}
		}
		row++;
	}

	contents_.clear();
	scalar_contents_.clear();
	z_ = new NCPA::VectorWithUnits( nlines, columns[ 0 ], depunits );
	
	for (i = 0; i < keys.size(); i++) {
		if (column_numbers[ i ] == 0) {
			this->add_property( keys[ i ], values[ i ], units[ i ] );
		} else if (column_numbers[ i ] > 1) {
			this->add_property( keys[ i ], nlines, columns[ column_numbers[ i ] - 1 ], units[ i ] );
		}
	}

	for ( i = 0; i < ncols; i++ ) {
		delete [] columns[ i ];
	}
}


NCPA::Atmosphere1D::~Atmosphere1D() {
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it=contents_.begin(); it != contents_.end(); ++it) {
		delete it->second;
	}
	for (std::map< std::string, NCPA::ScalarWithUnits * >::iterator it=scalar_contents_.begin(); 
			it != scalar_contents_.end(); ++it) {
		delete it->second;
	}
	contents_.clear();
	scalar_contents_.clear();
	delete z_;
}

void NCPA::Atmosphere1D::read_attenuation_from_file( const std::string &new_key,
	const std::string &filename ) {

	std::ifstream in( filename );
	std::string line;
	std::ostringstream oss;
	std::getline( in, line );
	std::vector< std::string > atmlines, fields;
	size_t i;

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
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
			oss << "Atmosphere1D - Error parsing attenuation line:" << std::endl << line << std::endl 
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

	// extend if necessary
	NCPA::units_t alt_units = get_altitude_units();
	this->convert_altitude_units( NCPA::Units::fromString( "km" ) );
	bool short_on_bottom = (z_a[0] > get_minimum_altitude());
	bool short_on_top    = (z_a[nlines-1] < get_maximum_altitude());
	double *temp_a, *temp_z;
	if (short_on_bottom) {
		temp_z = new double[ nlines+1 ];
		temp_a = new double[ nlines+1 ];
		std::memcpy( temp_z+1, z_a, nlines );
		std::memcpy( temp_a+1, attn, nlines );
		temp_z[ 0 ] = get_minimum_altitude();
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
		temp_z[ nlines-1 ] = get_maximum_altitude();
		temp_a[ nlines-1 ] = temp_a[ nlines-2 ];
		delete [] z_a;
		delete [] attn;
		z_a = temp_z;
		attn = temp_a;
		nlines++;
	}

	// Now interpolate onto current z grid
	double *existing_z = new double[ this->nz() ];
	this->get_altitude_vector( existing_z );
	double *interp_attn = new double[ this->nz() ];
	std::memset( interp_attn, 0, this->nz() * sizeof( double ) );

	gsl_interp_accel *aacc = gsl_interp_accel_alloc();
	gsl_spline *aspl = gsl_spline_alloc( gsl_interp_cspline, nlines );
	gsl_spline_init( aspl, z_a, attn, nlines );
	for (i = 0; i < this->nz(); i++) {
		interp_attn[ i ] = gsl_spline_eval( aspl, existing_z[ i ], aacc );
	}
	gsl_spline_free( aspl );
	gsl_interp_accel_free( aacc );

	this->add_property( new_key, this->nz(), interp_attn, NCPA::Units::fromString( "km" ) );
	this->convert_altitude_units( alt_units );
	delete [] interp_attn;
	delete [] existing_z;
	delete [] attn;
	delete [] z_a;

}


size_t NCPA::Atmosphere1D::get_basis_length() const {
	//return nz_;
	return z_->size();
}

size_t NCPA::Atmosphere1D::nz() const {
	return get_basis_length();
}

void NCPA::Atmosphere1D::copy_vector_property( const std::string &old_key,
		const std::string &new_key ) {
	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}
	if (!contains_vector( old_key )) {
		throw std::runtime_error( "Requested source key " + old_key + " not found" );
	}

	NCPA::AtmosphericProperty1D *source = contents_.at( old_key );
	units_t source_units;
	double *source_data = new double[ nz() ];
	source->get_vector( source_data, &source_units );
	add_property( new_key, nz(), source_data, source_units );
	delete [] source_data;
}

NCPA::AtmosphericProperty1D *NCPA::Atmosphere1D::get_vector_property_object(
			const std::string &key ) const {
	if (!contains_vector( key )) {
		throw std::runtime_error( "Requested key " + key + " not found." );
	}
	return contents_.at( key );
}

NCPA::ScalarWithUnits *NCPA::Atmosphere1D::get_scalar_property_object(
			const std::string &key ) const {
	if (!contains_scalar( key )) {
		throw std::runtime_error( "Requested key " + key + " not found." );
	}
	return scalar_contents_.at( key );
}

void NCPA::Atmosphere1D::copy_scalar_property( const std::string &old_key,
		const std::string &new_key ) {
	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}
	if (!contains_scalar( old_key )) {
		throw std::runtime_error( "Requested source key " + old_key + " not found" );
	}

	NCPA::ScalarWithUnits *source = scalar_contents_.at( old_key );
	units_t source_units = source->get_units();
	double value = source->get();
	add_property( new_key, value, source_units );
}

void NCPA::Atmosphere1D::calculate_wind_speed( const std::string &new_key,
		const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}


	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *u_prop = contents_.at( we_wind_speed_key );
	NCPA::AtmosphericProperty1D *v_prop = contents_.at( sn_wind_speed_key );
	units_t u_units, v_units;
	double *u = new double[ nz_ ];
	double *v = new double[ nz_ ];
	u_prop->get_vector( u, &u_units );
	v_prop->get_vector( v, &v_units );
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed properties" );
	}

	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( u[ i ] * u[ i ] + v[ i ] * v[ i ] );
	}
	delete [] u;
	delete [] v;
	add_property( new_key, nz_, c, u_units );
	delete [] c;
}

void NCPA::Atmosphere1D::calculate_wind_direction( const std::string &new_key,
	const std::string &we_wind_speed_key,
	const std::string &sn_wind_speed_key, units_t direction_units ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}

	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *u_prop = contents_.at( we_wind_speed_key );
	NCPA::AtmosphericProperty1D *v_prop = contents_.at( sn_wind_speed_key );
	units_t u_units, v_units;
	double *u = new double[ nz_ ];
	double *v = new double[ nz_ ];
	u_prop->get_vector( u, &u_units );
	v_prop->get_vector( v, &v_units );
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed properties" );
	}

	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = PI/2.0 - std::atan2( v[ i ], u[ i ] );
	}
	delete [] u;
	delete [] v;
	NCPA::Units::convert( c, nz_, NCPA::UNITS_ANGLE_RADIANS, NCPA::UNITS_ANGLE_DEGREES, c );
	NCPA::Units::convert( c, nz_, NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, direction_units, c );
	add_property( new_key, nz_, c, direction_units );
	delete [] c;
}

void NCPA::Atmosphere1D::calculate_wind_component( const std::string &new_key,
	const std::string &wind_speed_key, const std::string &wind_direction_key,
	double azimuth ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}
	
	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *sp_prop = contents_.at( wind_speed_key );
	NCPA::AtmosphericProperty1D *dir_prop = contents_.at( wind_direction_key );
	//NCPA::ScalarWithUnits *prop_dir = scalar_contents_.at( prop_direction_key );
	units_t sp_units, dir_units;
	double *sp = new double[ nz_ ];
	double *dir = new double[ nz_ ];
	sp_prop->get_vector( sp, &sp_units );
	dir_prop->get_vector( dir, &dir_units );
	//double propaz_rad = NCPA::Units::convert( prop_dir->get(), prop_dir->get_units(), NCPA::UNITS_ANGLE_RADIANS );

	// convert wind direction to radians for trig
	NCPA::Units::convert( dir, nz_, NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS, dir );
	double az_rad = NCPA::Units::convert( azimuth, NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS );

	double *wc = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		wc[ i ] = sp[ i ] * std::cos( dir[ i ] - az_rad );
	}
	delete [] sp;
	delete [] dir;
	add_property( new_key, nz_, wc, sp_units );
	delete [] wc;
}

void NCPA::Atmosphere1D::calculate_effective_sound_speed( const std::string &new_key,
	const std::string &sound_speed_key, const std::string &wind_component_key ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}


	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *c0_prop = contents_.at( sound_speed_key );
	NCPA::AtmosphericProperty1D *ws_prop = contents_.at( wind_component_key );
	units_t c0_units, ws_units;
	double *c0 = new double[ nz_ ];
	double *ws = new double[ nz_ ];
	c0_prop->get_vector( c0, &c0_units );
	ws_prop->get_vector( ws, &ws_units );

	// make sure units are consistent
	if (c0_units != ws_units) {
		NCPA::Units::convert( ws, nz_, ws_units, c0_units, ws );
		ws_units = c0_units;
	}

	double *ceff = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		ceff[ i ] = ws[ i ] + c0[ i ];
	}
	delete [] ws;
	delete [] c0;
	add_property( new_key, nz_, ceff, c0_units );
	delete [] ceff;
}

void NCPA::Atmosphere1D::get_altitude_vector( double *buffer, units_t *buffer_units ) const {
	z_->get_vector( buffer, buffer_units );
}

NCPA::units_t NCPA::Atmosphere1D::get_altitude_units() const {
	return z_->get_units();
}

void NCPA::Atmosphere1D::get_altitude_vector( double *buffer ) const {
	z_->get_vector( buffer );
}

void NCPA::Atmosphere1D::get_property_vector( const std::string &key,
		double *buffer, units_t *buffer_units ) const {
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	prop->get_vector( buffer, buffer_units );
}

void NCPA::Atmosphere1D::get_property_vector( const std::string &key,
		double *buffer ) const {
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	prop->get_vector( buffer );
}

NCPA::units_t NCPA::Atmosphere1D::get_property_units( const std::string &key ) const {

	if (contents_.count( key ) == 1) {
		return contents_.at( key )->get_units();
	}
	if (scalar_contents_.count( key ) == 1) {
		return scalar_contents_.at( key )->get_units();
	}
	throw std::out_of_range( "No vector or scalar quantity found with key " + key );
}


double NCPA::Atmosphere1D::get_minimum_altitude() const {
	return (*z_)[ 0 ];
}

double NCPA::Atmosphere1D::get_maximum_altitude() const {
	return (*z_)[ get_basis_length() - 1 ];
}

void NCPA::Atmosphere1D::add_property( const std::string &key, size_t n_points,
		double *quantity_points, units_t quantity_units ) {

	// see if we already have one with that key
	if (contains_key( key )) {
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	}

	size_t nz_ = get_basis_length();
	if (n_points != nz_) {
		std::ostringstream oss;
		oss << "Array length " << n_points << " for key " << key << " does not match number of altitude points " << nz_;
		throw std::runtime_error( oss.str() );
	}
	double *z_values = new double[ nz_ ];
	NCPA::units_t z_units_ = NCPA::UNITS_NONE;
	z_->get_vector( z_values, &z_units_ );

	NCPA::AtmosphericProperty1D *prop = new AtmosphericProperty1D( n_points, z_values, z_units_, quantity_points, quantity_units );
	contents_[ key ] = prop;
	delete [] z_values;
}

void NCPA::Atmosphere1D::add_property( const std::string &key, double value,
		NCPA::units_t units ) {
	if (contains_key( key )) {
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	}
	NCPA::ScalarWithUnits *scalar = new ScalarWithUnits( value, units );
	scalar_contents_[ key ] = scalar;
}

double NCPA::Atmosphere1D::get( const std::string &key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	if (contains_vector( key )) {
		return contents_.at( key )->get( altitude );
	} else {
		throw std::out_of_range( "No vector quantity \"" + key + "\" found" );
	}
}

double NCPA::Atmosphere1D::get( const std::string &key ) const {

	if (contains_scalar( key )) {
		return scalar_contents_.at( key )->get();
	} else {
		throw std::out_of_range( "No scalar quantity \"" + key + "\" found" );
	}

}

bool NCPA::Atmosphere1D::contains_vector( const std::string &key ) const {
	return (contents_.count( key ) == 1);
}

bool NCPA::Atmosphere1D::contains_scalar( const std::string &key ) const {
	return (scalar_contents_.count( key ) == 1);
}

bool NCPA::Atmosphere1D::contains_key( const std::string &key ) const {
	return contains_vector( key ) || contains_scalar( key );
}

double NCPA::Atmosphere1D::get_first_derivative( const std::string &key,
		double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_first_derivative( altitude );
}


double NCPA::Atmosphere1D::get_second_derivative( const std::string &key,
		double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_second_derivative( altitude );
}


void NCPA::Atmosphere1D::calculate_sound_speed_from_temperature(
	const std::string &new_key, const std::string &temperature_key,
	NCPA::units_t wind_units ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}

	NCPA::AtmosphericProperty1D *t_prop = contents_.at( temperature_key );
	NCPA::units_t old_units = t_prop->get_units();
	t_prop->convert_units( NCPA::UNITS_TEMPERATURE_KELVIN );
	size_t nz_ = get_basis_length();
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		//c[ i ] = std::sqrt( GAMMA_FOR_C * R_FOR_C * t_prop->get( (*z_)[i] ) );
		c[ i ] = NCPA::AtmosphericModel::soundspeed_from_temperature(
			t_prop->get( (*z_)[i] ) );
	}
	t_prop->convert_units( old_units );

	NCPA::Units::convert( c, nz_, NCPA::UNITS_SPEED_METERS_PER_SECOND, wind_units, c );
	add_property( new_key, nz_, c, wind_units );
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_pressure_and_density(
		const std::string &new_key, const std::string &pressure_key,
		const std::string &density_key, NCPA::units_t wind_units ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}

	NCPA::AtmosphericProperty1D *p_prop = contents_.at( pressure_key );
	NCPA::units_t old_p_units = p_prop->get_units();
	p_prop->convert_units( NCPA::UNITS_PRESSURE_PASCALS );
	NCPA::AtmosphericProperty1D *r_prop = contents_.at( density_key );
	NCPA::units_t old_r_units = r_prop->get_units();
	r_prop->convert_units( NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER );
	size_t nz_ = get_basis_length();
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = NCPA::AtmosphericModel::soundspeed_from_pressure_density(
			p_prop->get( (*z_)[i] ), r_prop->get( (*z_)[i] ) );
		// c[ i ] = std::sqrt( GAMMA_FOR_C * p_prop->get( (*z_)[i] ) / r_prop->get( (*z_)[i] ) );
	}

	p_prop->convert_units( old_p_units );
	r_prop->convert_units( old_r_units );
	NCPA::Units::convert( c, nz_, NCPA::UNITS_SPEED_METERS_PER_SECOND, wind_units, c );
	add_property( new_key, nz_, c, wind_units );
	delete [] c;
}

void NCPA::Atmosphere1D::convert_altitude_units( units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	z_->convert_units( new_units );

	// update all contents
	for ( std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it = contents_.begin();
						it != contents_.end(); ++it ) {
		(*it).second->convert_altitude_units( new_units );
	}
}

void NCPA::Atmosphere1D::convert_property_units( const std::string &key,
		units_t new_units ) {
	// if it's not there it'll throw an out_of_range exception
	if (contains_vector( key )) {
		contents_.at( key )->convert_units( new_units );
		return;
	}
	if (contains_scalar( key )) {
		scalar_contents_.at( key )->convert_units( new_units );
		return;
	}
	throw std::out_of_range( "No vector or scalar quantity found with key " + key );	
}



void NCPA::Atmosphere1D::do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double *units_buffer = new double[ n_points ];
	std::memset( units_buffer, 0, n_points * sizeof( double ) );
	
	// throws out_of_range if conversion is undefined
	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );

	// successful, so record the units change
	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
	delete [] units_buffer;
}

std::vector< std::string > NCPA::Atmosphere1D::get_vector_keys() const {
	std::vector< std::string > keys;
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::const_iterator it = contents_.begin();
			it != contents_.end(); ++it) {
		keys.push_back( (*it).first );
	}
	return keys;
}

std::vector< std::string > NCPA::Atmosphere1D::get_scalar_keys() const {
	std::vector< std::string > keys;
	for (std::map< std::string, NCPA::ScalarWithUnits * >::const_iterator it = scalar_contents_.begin();
			it != scalar_contents_.end(); ++it) {
		keys.push_back( (*it).first );
	}
	return keys;
}

std::vector< std::string > NCPA::Atmosphere1D::get_keys() const {
	std::vector< std::string > vkeys = get_vector_keys();
	std::vector< std::string > skeys = get_scalar_keys();
	std::vector< std::string > allkeys( vkeys );
	allkeys.reserve( vkeys.size() + skeys.size() );
	for (std::vector<std::string>::const_iterator cit = skeys.begin();
		cit != skeys.end(); ++cit) {
		allkeys.push_back( *cit );
	}
	return allkeys;
}

void NCPA::Atmosphere1D::print_atmosphere( const std::vector< std::string >& columnorder,
	const std::string &altitude_key, std::ostream& os ) {

	// check columnorder variable for key validity
	std::vector< std::string >::const_iterator vit;
	for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit) {
		if (! contains_vector( *vit ) ) {
			throw std::invalid_argument( "No vector quantity exists with key " + *vit );
		}
	}

	// first we do the header.  That contains column descriptions as well as scalar values
	// scalars first
	for (auto mit = scalar_contents_.cbegin(); mit != scalar_contents_.cend(); ++mit ) {
		os  << "#% 0, " << (*mit).first << ", " 
			<< NCPA::Units::toStr( (*mit).second->get_units() ) << ", "
			<< (*mit).second->get() << std::endl;
	}

	// Now column descriptors.  Altitude first
	os  << "#% 1, " << altitude_key << ", " 
		<< NCPA::Units::toStr( z_->get_units() ) << std::endl;
	unsigned int column = 2;
	for ( vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
		os  << "#% " << column << ", "
			<< *vit << ", "
			<< NCPA::Units::toStr( get_property_units( *vit ) ) << std::endl;
		column++;
	}

	// Now columns
	size_t nz_ = get_basis_length();
	os.setf( std::ios::scientific, 	std::ios::floatfield );
	os.setf( std::ios::right, 		std::ios::adjustfield );
	os.precision( 6 );
	os.width( 9 );
	os.fill( ' ' );
	for ( size_t i = 0; i < nz_; i++) {
		os << (*z_)[ i ];
		for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
			os << " " << get( *vit, (*z_)[ i ] );
		}
		os << std::endl;
	}
	os.flush();
}

void NCPA::Atmosphere1D::print_atmosphere( const std::string &altitude_key,
		std::ostream& os ) {
	print_atmosphere( get_keys(), altitude_key, os );
}

void NCPA::Atmosphere1D::resample( double new_dz ) {

	// resample altitude
	size_t old_nz = get_basis_length();
	double *old_z = new double[ old_nz ];
	NCPA::units_t old_z_units;
	z_->get_vector( old_z, &old_z_units );
	double z0 = old_z[0];
	double z1 = old_z[ old_nz - 1 ];
	size_t new_nz = (size_t)std::floor( ( z1 - z0 ) / new_dz ) + 1;

	double *new_z = new double[ new_nz ];
	//double *new_prop = new double[ new_nz ];
	for (size_t i = 0; i < new_nz; i++) {
		new_z[ i ] = z0 + ((double)i) * new_dz;
		//new_prop[ i ] = get( new_z[ i ] );
	}

	delete z_;
	z_ = new NCPA::VectorWithUnits( new_nz, new_z, old_z_units );

	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it = contents_.begin();
			it != contents_.end(); ++it) {
		(*it).second->resample( new_dz );
	}

	delete [] new_z;
	delete [] old_z;
}

void NCPA::Atmosphere1D::remove_property( const std::string &key ) {

	if (contains_vector( key )) {
		NCPA::AtmosphericProperty1D *prop = contents_.at( key );
		delete prop;
		contents_.erase( key );
	}
	if (contains_scalar( key )) {
		NCPA::ScalarWithUnits *swu = scalar_contents_.at( key );
		delete swu;
		scalar_contents_.erase( key );
	}
}


void NCPA::Atmosphere1D::calculate_attenuation( const std::string &new_key,
		const std::string &temperature_key, const std::string &pressure_key,
		const std::string &density_key, double freq, double tweak ) {

	if (contains_key( new_key )) {
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	}
	
	size_t nz = this->nz();
	double *Z = new double[ nz ];
	double *T = new double[ nz ];
	double *P = new double[ nz ];
	double *D = new double[ nz ];
	double *A = new double[ nz ];
	NCPA::units_t old_z_units, old_t_units, old_p_units, old_d_units;
	get_altitude_vector( Z, &old_z_units );
	get_property_vector( temperature_key, 	T, &old_t_units );
	get_property_vector( pressure_key, 		P, &old_p_units );
	get_property_vector( density_key, 		D, &old_d_units );

	// we do the units conversion in-place on the raw vectors so we don't have to change them back later
	NCPA::Units::convert( Z, nz, old_z_units, NCPA::Units::fromString( "km" ), Z );
	NCPA::Units::convert( T, nz, old_t_units, NCPA::Units::fromString( "K" ), T );
	NCPA::Units::convert( P, nz, old_p_units, NCPA::Units::fromString( "Pa" ), P );
	NCPA::Units::convert( D, nz, old_d_units, NCPA::Units::fromString( "kg/m3" ), D );

	// set up constants for attenuation calculation
	// double mu_o  		= 18.192E-6;    		// Reference viscosity [kg/(m*s)]
	// // double T_o   		= T[0];         		// Reference temperature [K]
	// double T_o   		= 293.15;         		// Reference temperature [K]
	// // double P_o   		= P[0];        			// Reference pressure [Pa]
	// double P_o   		= 101325;        		// Reference pressure [Pa]
	// double S     		= 117.0;	     		// Sutherland constant [K]

	// // heat capacity|volume for O2, N2, CO2, and O3
	// double Cv_R[] = { 5.0/2.0, 5.0/2.0, 3.0, 3.0 };

	// // heat capacity|pressure for O2, N2, CO2, and O3
	// double Cp_R[] = { 7.0/2.0, 7.0/2.0, 4.0, 4.0 };

	// // characteristic temperature for O2, N2, CO2, and O3
	// double theta[] = { 2239.1, 3352.0, 915.0, 1037.0 };

	for ( size_t ii = 0; ii < nz; ii++ ) {
		A[ii] = NCPA::AtmosphericModel::attenuation_sutherland_bass(
			Z[ii], T[ii], P[ii], D[ii], freq ) * tweak;
	}

	add_property( new_key, nz, A, NCPA::Units::fromString( "" ) );

	delete [] Z;
	delete [] T;
	delete [] P;
	delete [] D;
	delete [] A;
		// double z       = Z[ii];	// km	AGL
		// double T_z     = T[ii];	// K
		// double P_z     = P[ii];	// Pa;
		// double D_z		= D[ii];	// kg/m3
		// double c_snd_z = std::sqrt( GAMMA_FOR_C * P_z / D_z );  // m/s
		// double mu      = mu_o * std::sqrt( T_z/T_o ) * ( (1.0 + S/T_o) / (1.0 + S/T_z) ); // Viscosity [kg/(m*s)]
		// double nu      = ( 8.0 * PI * freq * mu ) / ( 3 * P_z );                   // Nondimensional frequency

		// // Gas fractions
		// double X[7];

		// //-------- Gas fraction polynomial fits -----------------------------------
		// if (z > 90.0) {                                         // O2 profile
		// 	X[0] = std::pow( 10,
		// 				49.296 - (1.5524*z) + (1.8714E-2*std::pow(z,2))
		// 		   		- (1.1069E-4*std::pow(z,3)) + (3.199E-7*std::pow(z,4))
		// 		   		- (3.6211E-10*std::pow(z,5))
		// 		);
		// } else {
		// 	X[0] = std::pow(10,-0.67887);
		// }

		// if (z > 76.0) {                                        // N2 profile
		// 	X[1] = std::pow(10,(1.3972E-1)-(5.6269E-3*z) + (3.9407E-5*std::pow(z,2) )
		// 		   - (1.0737E-7*std::pow(z,3)));
		// } else {
		// 	X[1] = std::pow(10,-0.10744);
		// }

		// X[2] = std::pow(10,-3.3979);                              // CO2 profile

		// if (z > 80.0) {                                        // O3 profile
		// 	X[3] = std::pow(10,-4.234-(3.0975E-2*z));
		// } else {
		// 	X[3] = std::pow(10,-19.027+(1.3093*z) - (4.6496E-2*std::pow(z,2))
		// 		   + (7.8543E-4*std::pow(z,3)) - (6.5169E-6*std::pow(z,4))
		// 		   + (2.1343E-8*std::pow(z,5)));
		// }

		// if (z > 95.0 ) {                                       // O profile
		// 	X[4] = std::pow(10,-3.2456+(4.6642E-2*z)-(2.6894E-4*std::pow(z,2))+(5.264E-7*std::pow(z,3)));
		// } else {
		// 	X[4] = std::pow(10,-11.195+(1.5408E-1*z)-(1.4348E-3*std::pow(z,2))+(1.0166E-5*std::pow(z,3)));
		// }

		// // N profile
		// X[5]  = std::pow(10,-53.746+(1.5439*z)-(1.8824E-2*std::pow(z,2))+(1.1587E-4*std::pow(z,3))
		// 	    -(3.5399E-7*std::pow(z,4))+(4.2609E-10*std::pow(z,5)));

		// if (z > 30.0 ) {                                        // H2O profile
		// 	X[6] = std::pow(10,-4.2563+(7.6245E-2*z)-(2.1824E-3*std::pow(z,2))-(2.3010E-6*std::pow(z,3))
		// 		   +(2.4265E-7*std::pow(z,4))-(1.2500E-09*std::pow(z,5)));
		// } else {
		// 	if (z > 100.) {
		// 		X[6] = std::pow(10,-0.62534-(8.3665E-2*z));
		// 	} else {
		// 		X[6] = std::pow(10,-1.7491+(4.4986E-2*z)-(6.8549E-2*std::pow(z,2))
		// 			   +(5.4639E-3*std::pow(z,3))-(1.5539E-4*std::pow(z,4))
		// 			   +(1.5063E-06*std::pow(z,5)));
		// 	}
		// }
		// double X_ON = (X[0] + X[1])/0.9903;

		// //-------- Rotational collision number-------------------------------------
		// double Z_rot_0 = 54.1*std::exp(-17.3*(std::pow(T_z,-1.0/3.0)));   // O2
		// double Z_rot_1 = 63.3*std::exp(-16.7*(std::pow(T_z,-1.0/3.0)));   // N2
		// double Z_rot_   = 1.0 / ( (X[1] / Z_rot_1) + (X[0] / Z_rot_0) );

		// //-------- Nondimensional atmospheric quantities---------------------------
		// double sigma = 5.0 / std::sqrt(21.0);
		// double nn    = (4.0/5.0) * std::sqrt(3.0/7.0) * Z_rot_;
		// double chi   = 3.0 * nn * nu / 4.0;
		// double cchi  = 2.36 * chi;

		// //---------Classical + rotational loss/dispersion--------------------------
		// double a_cl    = (2 * PI * freq / c_snd_z)
		// 					* std::sqrt( 0.5 * (std::sqrt(1+std::pow(nu,2))-1) * (1+std::pow(cchi,2))
		// 					/ ((1+std::pow(nu,2))*(1+std::pow(sigma*cchi,2))) );
		// double a_rot   = (2*PI*freq/c_snd_z) * X_ON
		// 					* ((pow(sigma,2)-1)*chi/(2*sigma))
		// 					* sqrt(0.5*(sqrt(1+pow(nu,2))+1)/((1+pow(nu,2))*(1+pow(cchi,2))));
		// double a_diff  = 0.003*a_cl;

		// //---------Vibrational relaxation-------------------------------------------
		// double Tr = std::pow( T_z/T_o ,-1.0/3.0 ) - 1.0;
		// double A1 = (X[0]+X[1]) * 24.0 * std::exp(-9.16*Tr);
		// double A2 = (X[4]+X[5]) * 2400.0;
		// double B  = 40400.0 * std::exp( 10.0*Tr );
		// double C  = 0.02 * std::exp( -11.2*Tr );
		// double D  = 0.391 * std::exp( 8.41*Tr );
		// double E  = 9.0 * std::exp( -19.9*Tr );
		// double F  = 60000.0;
		// double G  = 28000.0 * std::exp( -4.17*Tr );
		// double H  = 22000.0 * std::exp( -7.68*Tr );
		// double I  = 15100.0 * std::exp( -10.4*Tr );
		// double J  = 11500.0 * std::exp( -9.17*Tr );
		// double K  = (8.48E08) * std::exp( 9.17*Tr );
		// double L  = std::exp( -7.72*Tr );
		// double ZZ = H*X[2] + I*(X[0] + 0.5*X[4]) + J*(X[1] + 0.5*X[5] ) + K*(X[6] + X[3]);
		// double hu = 100.0 * (X[3] + X[6]);
		// double f_vib[ 4 ], a_vib_c[ 4 ];
		// f_vib[0] = (P_z/P_o) * (mu_o/mu) * (A1 + A2 + B*hu*(C+hu)*(D+hu));
		// f_vib[1] = (P_z/P_o) * (mu_o/mu) * (E + F*X[3] + G*X[6]);
		// f_vib[2] = (P_z/P_o) * (mu_o/mu) * ZZ;
		// f_vib[3] = (P_z/P_o) * (mu_o/mu) * (1.2E5)*L;

		// double a_vib = 0.0;
		// for (size_t m=0; m<4; m++)
		// {
		// 	double C_R        = ((std::pow(theta[m]/T_z,2))*std::exp(-theta[m]/T_z))/(std::pow(1-std::exp(-theta[m]/T_z),2));
		// 	double A_max      = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
		// 	       A_max      = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
		// 	a_vib_c[m] = (A_max/c_snd_z)*((2*(pow(freq,2))/f_vib[m])/(1+pow(freq/f_vib[m],2)));
		// 	a_vib      = a_vib + a_vib_c[m];
		// }

		// A[ii] = (a_cl + a_rot + a_diff + a_vib) * tweak;
	// }

	// add_property( new_key, nz, A, NCPA::Units::fromString( "" ) );

	// delete [] Z;
	// delete [] T;
	// delete [] P;
	// delete [] D;
	// delete [] A;
}
