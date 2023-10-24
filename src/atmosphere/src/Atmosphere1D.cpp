#include "Atmosphere1D.h"
#include "AtmosphericProperty1D.h"
#include "AtmosphericModel.h"
#include "NCPAUnits.h"
#include "NCPACommon.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "gsl/gsl_version.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

// #ifndef GAMMA_FOR_C
// #define GAMMA_FOR_C 1.4
// #endif

// #ifndef R_FOR_C
// #define R_FOR_C 287.0
// #endif

NCPA::Atmosphere1D::Atmosphere1D() : NCPA::AtmosphericModel() {}

NCPA::Atmosphere1D::Atmosphere1D( size_t n_altitude_points, double *altitude_points,
		units_t altitude_units ) : NCPA::AtmosphericModel() {
	reset_altitude_vector( n_altitude_points, altitude_points, altitude_units );
}

NCPA::Atmosphere1D::Atmosphere1D( size_t n_altitude_points, double *altitude_points,
		const std::string &altitude_units )
		: NCPA::Atmosphere1D::Atmosphere1D(
				n_altitude_points, altitude_points, NCPA::Units::fromString(altitude_units)) {}

NCPA::Atmosphere1D::Atmosphere1D( const std::string &filename,
	const std::string &headerfilename ) : NCPA::AtmosphericModel() {

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

	headerlines_.clear();
}

NCPA::Atmosphere1D::Atmosphere1D( const Atmosphere1D &source ) : NCPA::AtmosphericModel() {
	z_ = NCPA::VectorWithUnits( source.z_ );
	for ( auto it = source.contents_.cbegin(); it != source.contents_.cend(); ++it ) {
		contents_[ it->first ] = new NCPA::AtmosphericProperty1D( *(it->second) );
	}
	for ( auto it = source.scalar_contents_.cbegin(); it != source.scalar_contents_.cend(); ++it ) {
		scalar_contents_[ it->first ] = new NCPA::ScalarWithUnits( *(it->second) );
	}
}

NCPA::Atmosphere1D::Atmosphere1D( Atmosphere1D &&source ) noexcept : NCPA::AtmosphericModel() {
	swap( *this, source );
}

NCPA::Atmosphere1D &NCPA::Atmosphere1D::operator=( NCPA::Atmosphere1D a ) {
	swap( *this, a );
	return *this;
}

void swap( NCPA::Atmosphere1D &a, NCPA::Atmosphere1D &b ) noexcept {
	swap(static_cast<NCPA::AtmosphericModel&>(a), static_cast<NCPA::AtmosphericModel&>(b) );
	swap( a.contents_, b.contents_ );
	swap( a.scalar_contents_, b.scalar_contents_ );
	swap( a.z_, b.z_ );
	swap( a.headerlines_, b.headerlines_ );
}


void NCPA::Atmosphere1D::reset_altitude_vector( size_t n_altitude_points,
			double *altitude_points, units_t altitude_units ) {
	this->reset_altitude_vector( NCPA::VectorWithUnits(
			n_altitude_points, altitude_points, altitude_units ) );
}

void NCPA::Atmosphere1D::reset_altitude_vector( size_t n_altitude_points,
			double *altitude_points, const std::string &altitude_units ) {
	this->reset_altitude_vector( NCPA::VectorWithUnits(
			n_altitude_points, altitude_points, altitude_units ) );
}

void NCPA::Atmosphere1D::reset_altitude_vector( const NCPA::VectorWithUnits &v ) {
	NCPA::VectorWithUnits new_z( v );
	// make sure it's distance
	NCPA::units_t oldu = new_z.get_units();
	if (!NCPA::Units::can_convert(oldu,"m")) {
		throw std::logic_error("New altitude vector is not in distance units!");
	}
	new_z.convert_units(oldu);
	this->cleanup_vectors_();
	z_ = new_z;
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
				headerlines_.push_back( line.substr( 2 ) );
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
			// 	headerlines_.push_back( line.substr( 2 ) );
			// } // otherwise it's a regular comment and can be ignored

		} else if (line.size() == 0) {
			// skip empty lines
		} else {
			atmlines.push_back( line );
		}

		std::getline( in, line );
	}
	//in.close();
	//cout << "Found " << headerlines_.size() << " header lines" << endl;
	//cout << "Found " << atmlines.size() << " data lines" << endl;

	// parse them out
	size_t nfields = headerlines_.size();
	if (nfields == 0) {
		throw std::runtime_error( "Atmosphere1D - No descriptive fields found." );
	}

	// hold contents
	std::vector< std::string > keys, fields;
	std::vector< unsigned int > column_numbers;
	std::vector< double > values;
	std::vector< NCPA::units_t > units;

	for (std::vector< std::string >::const_iterator it = headerlines_.begin(); it != headerlines_.end(); ++it) {
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
		//NCPA::units_t tempunits = NCPA::units_t::NONE;
		NCPA::units_t tempunits = NCPA::Units::fromString( NCPA::deblank( fields[ 2 ] ) );
		if (tempunits == NCPA::units_t::NONE) {
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
			try {
				tempval = std::stod( NCPA::deblank(fields[ 3 ]) );
			} catch (std::out_of_range &e) {
				oss << "Error converting data line:" << std::endl << *it
					<< std::endl << "Value " << fields[3]
					<< " is out of range for double precision.";
				throw std::runtime_error( oss.str() );
			}
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
	NCPA::units_t depunits = NCPA::units_t::NONE;
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
				thiscol[ row ] = std::stod( fields[ i ] );
			} catch (std::out_of_range &e) {
				oss << "Error converting data line:" << std::endl << *it
					<< std::endl << "Value " << fields[i]
					<< " is out of range for double precision.";
				for (i = 0; i < ncols; i++) {
					delete [] columns[i];
				}
				throw std::runtime_error( oss.str() );
			} catch (std::invalid_argument &e) {
				oss << "Error parsing data line:" << std::endl << *it << std::endl
					<< "Can't parse field " << fields[ i ] << " as a double";
				for ( i = 0; i < ncols; i++ ) {
					delete [] columns[ i ];
				}
				throw std::runtime_error( oss.str() );
			}
		}
		row++;
	}

//	contents_.clear();
//	scalar_contents_.clear();
	this->cleanup_();
	z_ = NCPA::VectorWithUnits( nlines, columns[ 0 ], depunits );
	
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
	cleanup_();
}

void NCPA::Atmosphere1D::cleanup_() {
	this->cleanup_vectors_();
	this->cleanup_scalars_();
	z_.clear();
}

void NCPA::Atmosphere1D::cleanup_vectors_() {
	for (auto it=contents_.begin(); it != contents_.end(); ++it) {
		delete it->second;
	}
	contents_.clear();
}

void NCPA::Atmosphere1D::cleanup_scalars_() {
	for (auto it=scalar_contents_.begin(); it != scalar_contents_.end(); ++it) {
		delete it->second;
	}
	scalar_contents_.clear();
}

void NCPA::Atmosphere1D::reset_splines() {
	for (auto it=contents_.begin(); it != contents_.end(); ++it) {
		it->second->reset_splines();
	}
}

void NCPA::Atmosphere1D::read_attenuation_from_file( const std::string &new_key,
	const std::string &filename ) {

	std::ifstream in( filename );
	std::string line;
	std::ostringstream oss;
	std::getline( in, line );
	std::vector< std::string > atmlines, fields;
	size_t i;

	assert_key_does_not_exist_(new_key);

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
			this_z = std::stod( fields[ 0 ] );
			this_a = std::stod( fields[ 1 ] );
		} catch ( std::invalid_argument &e ) {
			oss << "Atmosphere1D - Error parsing attenuation line:" << std::endl << line << std::endl 
				<< "Both fields must be numerical" << std::endl;
			throw std::invalid_argument( oss.str() );
		} catch (std::out_of_range &e) {
			oss << "Error converting data line:" << std::endl << line
				<< std::endl << "Value  is out of range for double precision.";
			throw std::runtime_error( oss.str() );
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
	gsl_spline *aspl = gsl_spline_alloc( ATMOSPHERIC_INTERPOLATION_TYPE, nlines );
//#if GSL_MAJOR_VERSION > 1
//	gsl_spline *aspl = gsl_spline_alloc( gsl_interp_steffen, nlines );
//#else
//	gsl_spline *aspl = gsl_spline_alloc( gsl_interp_cspline, nlines );
//#endif
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
	return z_.size();
}

size_t NCPA::Atmosphere1D::nz() const {
	return get_basis_length();
}

void NCPA::Atmosphere1D::copy_vector_property( const std::string &old_key,
		const std::string &new_key ) {
	assert_key_does_not_exist_(new_key);

	NCPA::AtmosphericProperty1D *source = get_vector_property_object( old_key );
	units_t source_units;
	double *source_data = new double[ nz() ];
	source->as_array( source_data, source_units );
	add_property( new_key, nz(), source_data, source_units );
	delete [] source_data;
}

NCPA::AtmosphericProperty1D *NCPA::Atmosphere1D::get_vector_property_object(
			const std::string &key ) {
	assert_vector_key_exists_(key);
	return contents_.at( key );
}

const NCPA::AtmosphericProperty1D * const NCPA::Atmosphere1D::get_const_vector_property_object(
			const std::string &key ) const {
	assert_vector_key_exists_(key);
	return contents_.at( key );
}

NCPA::ScalarWithUnits *NCPA::Atmosphere1D::get_scalar_property_object(
			const std::string &key ) {
	assert_scalar_key_exists_(key);
	return scalar_contents_.at( key );
}

const NCPA::ScalarWithUnits * const NCPA::Atmosphere1D::get_const_scalar_property_object(
			const std::string &key ) const {
	assert_scalar_key_exists_(key);
	return scalar_contents_.at( key );
}

void NCPA::Atmosphere1D::copy_scalar_property( const std::string &old_key,
		const std::string &new_key ) {
	assert_key_does_not_exist_(new_key);

	NCPA::ScalarWithUnits *source = this->get_scalar_property_object( old_key );
	units_t source_units = source->get_units();
	double value = source->get();
	add_property( new_key, value, source_units );
}

void NCPA::Atmosphere1D::calculate_wind_speed( const std::string &new_key,
		const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key ) {

	assert_key_does_not_exist_(new_key);


	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *u_prop = get_vector_property_object( we_wind_speed_key );
	NCPA::AtmosphericProperty1D *v_prop = get_vector_property_object( sn_wind_speed_key );
	units_t u_units, v_units;
	double *u = new double[ nz_ ];
	double *v = new double[ nz_ ];
	u_prop->as_array( u, u_units );
	v_prop->as_array( v, v_units );
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

	assert_key_does_not_exist_(new_key);

	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *u_prop = get_vector_property_object( we_wind_speed_key );
	NCPA::AtmosphericProperty1D *v_prop = get_vector_property_object( sn_wind_speed_key );
	units_t u_units, v_units;
	double *u = new double[ nz_ ];
	double *v = new double[ nz_ ];
	u_prop->as_array( u, u_units );
	v_prop->as_array( v, v_units );
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed properties" );
	}

	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
//		if (v[i] == 0.0 && u[i] == 0.0) {
//			c[ i ] = 0.0;
//		} else {
			c[ i ] = M_PI/2.0 - std::atan2( v[ i ], u[ i ] );
//		}
	}
	delete [] u;
	delete [] v;
	NCPA::Units::convert( c, nz_, NCPA::units_t::ANGLE_RADIANS, NCPA::units_t::ANGLE_DEGREES, c );
	NCPA::Units::convert( c, nz_, NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, direction_units, c );
	add_property( new_key, nz_, c, direction_units );
	delete [] c;
}

void NCPA::Atmosphere1D::calculate_wind_component( const std::string &new_key,
	const std::string &wind_speed_key, const std::string &wind_direction_key,
	double azimuth ) {

	assert_key_does_not_exist_(new_key);
	
	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *sp_prop = get_vector_property_object( wind_speed_key );
	NCPA::AtmosphericProperty1D *dir_prop = get_vector_property_object( wind_direction_key );
	//NCPA::ScalarWithUnits *prop_dir = scalar_contents_.at( prop_direction_key );
	units_t sp_units, dir_units;
	double *sp = new double[ nz_ ];
	double *dir = new double[ nz_ ];
	sp_prop->as_array( sp, sp_units );
	dir_prop->as_array( dir, dir_units );
	//double propaz_rad = NCPA::Units::convert( prop_dir->get(), prop_dir->get_units(), NCPA::units_t::ANGLE_RADIANS );

	// convert wind direction to radians for trig
	NCPA::Units::convert( dir, nz_, NCPA::units_t::ANGLE_DEGREES, NCPA::units_t::ANGLE_RADIANS, dir );
	double az_rad = NCPA::Units::convert( azimuth, NCPA::units_t::ANGLE_DEGREES, NCPA::units_t::ANGLE_RADIANS );

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

	assert_key_does_not_exist_(new_key);


	size_t nz_ = get_basis_length();
	NCPA::AtmosphericProperty1D *c0_prop = get_vector_property_object( sound_speed_key );
	NCPA::AtmosphericProperty1D *ws_prop = get_vector_property_object( wind_component_key );
	units_t c0_units, ws_units;
	double *c0 = new double[ nz_ ];
	double *ws = new double[ nz_ ];
	c0_prop->as_array( c0, c0_units );
	ws_prop->as_array( ws, ws_units );

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

void NCPA::Atmosphere1D::get_altitude_vector( double *buffer, units_t &buffer_units ) const {
	NCPA::VectorWithUnits ztemp(z_);
	ztemp.as_array( buffer, buffer_units );
}

NCPA::units_t NCPA::Atmosphere1D::get_altitude_units() const {
	NCPA::VectorWithUnits ztemp(z_);
	return ztemp.get_units();
}

void NCPA::Atmosphere1D::get_altitude_vector( double *buffer ) const {
	NCPA::VectorWithUnits ztemp(z_);
	ztemp.get_values( buffer );
}

void NCPA::Atmosphere1D::get_property_vector( const std::string &key,
		double *buffer, units_t &buffer_units ) const {
	NCPA::AtmosphericProperty1D vtemp(*(get_const_vector_property_object( key )));
	vtemp.get_vector_as( buffer_units ).as_array( buffer );
}

void NCPA::Atmosphere1D::get_property_vector( const std::string &key,
		double *buffer ) const {
	NCPA::AtmosphericProperty1D vtemp(*(get_const_vector_property_object( key )));
	vtemp.as_array( buffer );
}

NCPA::units_t NCPA::Atmosphere1D::get_property_units( const std::string &key ) const {

	if (contents_.count( key ) == 1) {
		return NCPA::AtmosphericProperty1D(
				*(get_const_vector_property_object( key ))
			).get_units();
	}
	if (scalar_contents_.count( key ) == 1) {
		return NCPA::ScalarWithUnits(
				*(get_const_scalar_property_object( key ))
			).get_units();
	}
	throw std::out_of_range( "No vector or scalar quantity found with key " + key );
}


double NCPA::Atmosphere1D::get_minimum_altitude() const {
	return z_.front().get();
}

double NCPA::Atmosphere1D::get_maximum_altitude() const {
	return z_.back().get();
}

void NCPA::Atmosphere1D::add_property( const std::string &key, size_t n_points,
		double *quantity_points, const std::string &quantity_units ) {
	this->add_property( key, NCPA::VectorWithUnits(
			n_points, quantity_points, quantity_units ) );
}

void NCPA::Atmosphere1D::add_property( const std::string &key, size_t n_points,
		double *quantity_points, units_t quantity_units ) {
	this->add_property( key,
			NCPA::VectorWithUnits(n_points, quantity_points, quantity_units) );
}

void NCPA::Atmosphere1D::add_property( const std::string &key,
		const NCPA::VectorWithUnits &v ) {
	assert_key_does_not_exist_(key);
	size_t nz_ = get_basis_length();
	if (v.size() != nz_) {
		std::ostringstream oss;
		oss << "Array length " << v.size() << " for key " << key
				<< " does not match number of altitude points " << nz_;
		throw std::runtime_error( oss.str() );
	}
	contents_[ key ] = new AtmosphericProperty1D( z_, v );
}

void NCPA::Atmosphere1D::add_property( const std::string &key, double value,
		NCPA::units_t units ) {
	this->add_property( key, NCPA::ScalarWithUnits( value, units ) );
}

void NCPA::Atmosphere1D::add_property( const std::string &key, double value,
		const std::string &units ) {
	this->add_property( key, NCPA::ScalarWithUnits( value, units ) );
}

void NCPA::Atmosphere1D::add_property( const std::string &key,
		const NCPA::ScalarWithUnits &s ) {
	assert_key_does_not_exist_(key);
	scalar_contents_[ key ] = new ScalarWithUnits( s );
}

double NCPA::Atmosphere1D::get( const std::string &key, double altitude ) const {
	// if it's not there it'll throw an out_of_range exception
	assert_vector_key_exists_(key);
	return get_const_vector_property_object( key )->get( altitude );
}

double NCPA::Atmosphere1D::get( const std::string &key ) const {
	assert_scalar_key_exists_(key);
	return get_const_scalar_property_object( key )->get();
}

double NCPA::Atmosphere1D::get_as( const std::string &key, NCPA::units_t asunits ) const {
	return this->get_const_scalar_property_object(key)->as(asunits);
}

double NCPA::Atmosphere1D::get_as( const std::string &key, double altitude, units_t u ) const {
	return this->get_const_vector_property_object(key)->get_as( altitude, u );
}

double NCPA::Atmosphere1D::get_first_derivative_as( const std::string &key, double altitude,
		units_t u ) const {
	return this->get_const_vector_property_object(key)->get_first_derivative_as( altitude, u );
}

double NCPA::Atmosphere1D::get_second_derivative_as( const std::string &key, double altitude,
		units_t u ) const {
	return this->get_const_vector_property_object(key)->get_second_derivative_as( altitude, u );
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

bool NCPA::Atmosphere1D::contains( const std::string &key ) const {
	return this->contains_key(key);
}

double NCPA::Atmosphere1D::get_first_derivative( const std::string &key,
		double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
//	NCPA::AtmosphericProperty1D *prop = get_vector_property_object( key );
	return get_const_vector_property_object( key )->get_first_derivative( altitude );
}


double NCPA::Atmosphere1D::get_second_derivative( const std::string &key,
		double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
//	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	return get_const_vector_property_object( key )->get_second_derivative( altitude );
}

void NCPA::Atmosphere1D::calculate_density_from_temperature_and_pressure(
	const std::string &new_key, const std::string &temperature_key,
	const std::string &pressure_key, const std::string &density_units ) {

	this->calculate_density_from_temperature_and_pressure(
			new_key, temperature_key, pressure_key, NCPA::Units::fromString(density_units));
}

void NCPA::Atmosphere1D::calculate_density_from_temperature_and_pressure(
	const std::string &new_key, const std::string &temperature_key,
	const std::string &pressure_key, NCPA::units_t density_units ) {

	assert_key_does_not_exist_(new_key);

	double *rho = new double[ z_.size() ];
	for (size_t i = 0; i < z_.size(); i++) {
		rho[ i ] = NCPA::AtmosphericModel::density_from_temperature_pressure(
			this->get( temperature_key, z_[i].get() ),
			this->get_property_units( temperature_key ),
			this->get( pressure_key, z_[i].get() ),
			this->get_property_units(pressure_key),
			density_units );
	}
	add_property( new_key, z_.size(), rho, density_units );
	delete [] rho;
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_temperature(
	const std::string &new_key, const std::string &temperature_key,
	const std::string &wind_units ) {
	this->calculate_sound_speed_from_temperature( new_key, temperature_key,
			NCPA::Units::fromString(wind_units) );
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_temperature(
	const std::string &new_key, const std::string &temperature_key,
	NCPA::units_t speed_units ) {

	assert_key_does_not_exist_(new_key);
	size_t nz_ = z_.size();
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = NCPA::AtmosphericModel::soundspeed_from_temperature(
				this->get( temperature_key, z_[i].get() ),
				this->get_property_units( temperature_key ),
				speed_units );
	}
	add_property( new_key, nz_, c, speed_units );
	delete [] c;
}

void NCPA::Atmosphere1D::calculate_temperature_from_sound_speed(
	const std::string &new_key, const std::string &speed_key,
	const std::string &temp_units ) {
	this->calculate_temperature_from_sound_speed(
			new_key,speed_key,NCPA::Units::fromString(temp_units));
}

void NCPA::Atmosphere1D::calculate_temperature_from_sound_speed(
	const std::string &new_key, const std::string &speed_key,
	NCPA::units_t temp_units ) {

	assert_key_does_not_exist_(new_key);
	size_t nz_ = z_.size();
	double *t = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		t[ i ] = NCPA::AtmosphericModel::temperature_from_soundspeed(
				this->get( speed_key, z_[i].get() ),
				this->get_property_units( speed_key ),
				temp_units );
	}
	add_property( new_key, nz_, t, temp_units );
	delete [] t;
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_pressure_and_density(
		const std::string &new_key, const std::string &pressure_key,
		const std::string &density_key, const std::string &wind_units ) {
	this->calculate_sound_speed_from_pressure_and_density(
			new_key, pressure_key, density_key, NCPA::Units::fromString(wind_units));
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_pressure_and_density(
		const std::string &new_key, const std::string &pressure_key,
		const std::string &density_key, NCPA::units_t speed_units ) {

	assert_key_does_not_exist_(new_key);
	size_t nz_ = z_.size();
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = NCPA::AtmosphericModel::soundspeed_from_pressure_density(
				this->get( pressure_key, z_[i].get() ),
				this->get_property_units( pressure_key ),
				this->get( density_key, z_[i].get() ),
				this->get_property_units( density_key ),
				speed_units );
	}
	add_property( new_key, nz_, c, speed_units );
	delete [] c;
}

void NCPA::Atmosphere1D::convert_altitude_units( units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	z_.convert_units( new_units );

	// update all contents
	for (auto it = contents_.begin(); it != contents_.end(); ++it) {
		it->second->convert_altitude_units( new_units );
	}
}

void NCPA::Atmosphere1D::convert_property_units( const std::string &key,
		units_t new_units ) {
	assert_key_exists_(key);
	// if it's not there it'll throw an out_of_range exception
	if (contains_vector( key )) {
		get_vector_property_object( key )->convert_units( new_units );
		return;
	}
	if (contains_scalar( key )) {
		get_scalar_property_object( key )->convert_units( new_units );
		return;
	}
//	throw std::out_of_range( "No vector or scalar quantity found with key " + key );
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
	for (auto it = contents_.begin(); it != contents_.end(); ++it) {
		keys.push_back( it->first );
	}
	return keys;
}

std::vector< std::string > NCPA::Atmosphere1D::get_scalar_keys() const {
	std::vector< std::string > keys;
	for (auto it = scalar_contents_.begin(); it != scalar_contents_.end(); ++it) {
		keys.push_back( it->first );
	}
	return keys;
}

std::vector< std::string > NCPA::Atmosphere1D::get_keys() const {
	std::vector< std::string > vkeys = get_vector_keys();
	std::vector< std::string > skeys = get_scalar_keys();
	std::sort(vkeys.begin(),vkeys.end());
	std::sort(skeys.begin(),skeys.end());
	std::vector< std::string > allkeys(vkeys.size() + skeys.size());
	std::vector< std::string >::iterator it = std::set_union(
			vkeys.begin(),vkeys.end(),skeys.begin(),skeys.end(),allkeys.begin());
	allkeys.resize(it-allkeys.begin());
	return allkeys;
}



void NCPA::Atmosphere1D::resample( double new_dz, NCPA::Interpolator1D *interp ) {

	// resample altitude
//	size_t old_nz = get_basis_length();
//	double *old_z = new double[ old_nz ];
//	NCPA::units_t old_z_units;
//	z_.as_array( old_z, old_z_units );
//	double z0 = old_z[0];
//	double z1 = old_z[ old_nz - 1 ];
//	size_t new_nz = (size_t)std::floor( ( z1 - z0 ) / new_dz ) + 1;

	double z0 = z_.front().get();
	double z1 = z_.back().get();
	size_t new_nz = (size_t)std::floor( ( z1 - z0 ) / new_dz ) + 1;
	double *new_z = new double[ new_nz ];
	for (size_t i = 0; i < new_nz; i++) {
		new_z[ i ] = z0 + ((double)i) * new_dz;
	}

	NCPA::VectorWithUnits nzv( new_nz, new_z, z_.get_units() );
	this->resample( nzv, interp );
}

void NCPA::Atmosphere1D::resample( NCPA::VectorWithUnits &new_z,
		NCPA::Interpolator1D *interp ) {

	NCPA::VectorWithUnits old_z( z_ );
	old_z.convert_units( new_z.get_units() );
	size_t 	old_nz = old_z.size(),
			new_nz = new_z.size();
	double *old_z_vals = NCPA::zeros<double>( old_nz );
	old_z.as_array( old_z_vals );

	// get copy of old vector map
	vector_atmospheric_property_map_t old_contents;
	swap( contents_, old_contents );
	this->reset_altitude_vector( new_z );  // resets contents at the same time

	interp->init();
	interp->allocate( old_nz );
	double *old_prop = NCPA::zeros<double>( old_nz );
	double *new_prop = NCPA::zeros<double>( new_nz );
	NCPA::units_t old_units;
	for (auto it = old_contents.begin(); it != old_contents.end(); ++it) {
		it->second->as_array( old_prop );
		interp->set( old_nz, old_z_vals, old_prop );
		for (size_t i = 0; i < new_nz; i++) {
			new_prop[i] = interp->f( new_z[i].get() );
		}
		NCPA::VectorWithUnits v( new_nz, new_prop, old_units );
		this->add_property( it->first, v );
		delete it->second;
		it->second = nullptr;
	}
	delete [] old_z_vals;
	delete [] old_prop;
	delete [] new_prop;
}

void NCPA::Atmosphere1D::remove_property( const std::string &key ) {

	if (contains_vector( key )) {
		NCPA::AtmosphericProperty1D *prop = get_vector_property_object( key );
		delete prop;
		contents_.erase( key );
	}
	if (contains_scalar( key )) {
		NCPA::ScalarWithUnits *swu = get_scalar_property_object( key );
		delete swu;
		scalar_contents_.erase( key );
	}
}


void NCPA::Atmosphere1D::calculate_attenuation( const std::string &new_key,
		const std::string &temperature_key, const std::string &pressure_key,
		const std::string &density_key, const std::string &humidity_key,
		double freq, double tweak ) {

	assert_key_does_not_exist_(new_key);

	// are we using density or humidity?
	if (density_key.size() > 0 && contains_vector( density_key )) {
		return calculate_attenuation( new_key, temperature_key,
			pressure_key, density_key, freq, tweak );
	} else if (humidity_key.size() > 0 && contains_vector( humidity_key )) {

		size_t nz = this->nz();
		double *Z = new double[ nz ];
		double *T = new double[ nz ];
		double *P = new double[ nz ];
		double *H = new double[ nz ];
		double *A = new double[ nz ];
		NCPA::units_t old_z_units, old_t_units, old_p_units;
		get_altitude_vector( Z, old_z_units );
		get_property_vector( temperature_key, 	T, old_t_units );
		get_property_vector( pressure_key, 		P, old_p_units );
		get_property_vector( humidity_key,      H );
		// get_property_vector( density_key, 		D, &old_d_units );

		// we do the units conversion in-place on the raw vectors so we don't have to change them back later
		NCPA::Units::convert( Z, nz, old_z_units, NCPA::Units::fromString( "km" ), Z );
		NCPA::Units::convert( T, nz, old_t_units, NCPA::Units::fromString( "K" ), T );
		NCPA::Units::convert( P, nz, old_p_units, NCPA::Units::fromString( "Pa" ), P );
		// NCPA::Units::convert( D, nz, old_d_units, NCPA::Units::fromString( "kg/m3" ), D );

		for ( size_t ii = 0; ii < nz; ii++ ) {
			A[ii] = NCPA::AtmosphericModel::attenuation_from_temperature_pressure_humidity(
				Z[ii], T[ii], P[ii], H[ii], freq ) * tweak;
		}

		add_property( new_key, nz, A, NCPA::Units::fromString( "" ) );

		delete [] Z;
		delete [] T;
		delete [] P;
		delete [] H;
		delete [] A;
	} else {
		throw std::runtime_error( "No valid density or humidity vector found!" );
	}
}


void NCPA::Atmosphere1D::calculate_attenuation( const std::string &new_key,
		const std::string &temperature_key, const std::string &pressure_key,
		const std::string &density_key, double freq, double tweak ) {

	assert_key_does_not_exist_(new_key);
	
	size_t nz = this->nz();
	double *Z = new double[ nz ];
	double *T = new double[ nz ];
	double *P = new double[ nz ];
	double *D = new double[ nz ];
	double *A = new double[ nz ];
	NCPA::units_t old_z_units, old_t_units, old_p_units, old_d_units;
	get_altitude_vector( Z, old_z_units );
	get_property_vector( temperature_key, 	T, old_t_units );
	get_property_vector( pressure_key, 		P, old_p_units );
	get_property_vector( density_key, 		D, old_d_units );

	// we do the units conversion in-place on the raw vectors so we don't have to change them back later
	NCPA::Units::convert( Z, nz, old_z_units, NCPA::Units::fromString( "km" ), Z );
	NCPA::Units::convert( T, nz, old_t_units, NCPA::Units::fromString( "K" ), T );
	NCPA::Units::convert( P, nz, old_p_units, NCPA::Units::fromString( "Pa" ), P );
	NCPA::Units::convert( D, nz, old_d_units, NCPA::Units::fromString( "kg/m3" ), D );

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
}

void NCPA::Atmosphere1D::assert_key_does_not_exist_( const std::string &key ) const {
	if (this->contains_key(key)) {
		throw std::invalid_argument( "Requested key " + key + " already exists in atmosphere" );
	}
}

void NCPA::Atmosphere1D::assert_key_exists_( const std::string &key ) const {
	if (!this->contains_key(key)) {
		throw std::out_of_range( "Requested key " + key + " does not exist in atmosphere" );
	}
}

void NCPA::Atmosphere1D::assert_vector_key_exists_( const std::string &key ) const {
	if (!this->contains_vector(key)) {
		throw std::out_of_range( "Requested key " + key + " does not exist in atmosphere" );
	}
}

void NCPA::Atmosphere1D::assert_scalar_key_exists_( const std::string &key ) const {
	if (!this->contains_scalar(key)) {
		throw std::out_of_range( "Requested key " + key + " does not exist in atmosphere" );
	}
}

void NCPA::Atmosphere1D::print_atmosphere(
				const std::vector< std::string >& properties,
				const std::string &altitude_key,
				std::ostream& os ) const {
	NCPA::VectorWithUnits tempz( z_ );
	this->print_atmosphere(tempz, properties, os, altitude_key );
}

void NCPA::Atmosphere1D::print_atmosphere(
		NCPA::VectorWithUnits &z_values,
		const std::vector< std::string >& properties,
		std::ostream& os, const std::string &z_label ) const {

	// check columnorder variable for key validity, and separate out
	// scalars from vectors
	std::vector< std::string > vectors, scalars;
	for (auto vit = properties.cbegin(); vit != properties.cend(); ++vit) {
		if (this->contains_vector( *vit )) {
			vectors.push_back( *vit );
		} else if (this->contains_scalar( *vit )) {
			scalars.push_back( *vit );
		} else {
			throw std::invalid_argument( "No quantity exists with key " + *vit );
		}
	}

	// Header.  Scalars first.
	for (auto sit = scalars.cbegin(); sit != scalars.cend(); ++sit) {
//		std::string u = NCPA::Units::toStr(this->get_property_units(*sit));
		os << this->make_scalar_header_line( *sit );
		os << std::endl;
	}

	// Now vectors.  Altitude first
	os << NCPA::Atmosphere1D::format_header_line( 1, z_label,
			NCPA::Units::toStr(z_values.get_units()) );
	os << std::endl;
	size_t column = 2;
	for (auto vit = vectors.cbegin(); vit != vectors.cend(); ++vit) {
		os << make_vector_header_line( *vit, column++ );
		os << std::endl;
	}

	// Now columns
	os.width( NCPA_ATMOSPHERE_FIELD_WIDTH );
	os.fill(' ');
	os.setf( std::ios::right, 		std::ios::adjustfield );
	for (auto zit = z_values.cbegin(); zit != z_values.cend(); ++zit) {
		os << NCPA::Atmosphere1D::format_for_stream(zit->get());
		for (auto vit = vectors.cbegin(); vit != vectors.cend(); ++vit) {
			os << " "
			   << NCPA::Atmosphere1D::format_for_stream(this->get( *vit, zit->get() ));
		}
		os << std::endl;
	}
//
//
//	os.setf( std::ios::scientific, 	std::ios::floatfield );
//	os.setf( std::ios::right, 		std::ios::adjustfield );
//	os.precision( 6 );
//	os.width( 9 );
//	os.fill( ' ' );
//	for ( size_t i = 0; i < nz_; i++) {
//		os << z_[ i ];
//		for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
//			os << " " << get( range, *vit, z_[ i ] );
//		}
//		os << std::endl;
//	}
	os.flush();
}

std::string NCPA::Atmosphere1D::make_scalar_header_line(
		const std::string &key ) const {
	std::ostringstream oss;
	if (!this->contains_scalar(key)) {
		oss << "No scalar property " << key << " in atmosphere!";
		throw std::invalid_argument( oss.str() );
	}
	return NCPA::Atmosphere1D::format_header_line(
			0,
			remove_underscores(key),
			NCPA::Units::toStr( this->get_property_units( key ) ),
			this->get(key) );
}

std::string NCPA::Atmosphere1D::make_vector_header_line(
		const std::string &key, size_t columnnumber ) const {
	std::ostringstream oss;
	if (!this->contains_vector(key)) {
		oss << "No vector property " << key << " in atmosphere!";
		throw std::invalid_argument( oss.str() );
	}
	return NCPA::Atmosphere1D::format_header_line(
			columnnumber, remove_underscores(key),
			NCPA::Units::toStr( this->get_property_units( key ) ) );
}

std::string NCPA::Atmosphere1D::make_header_line(
		const std::string &key, size_t columnnumber ) const {
	if (this->contains_vector(key)) {
		return this->make_vector_header_line( key, columnnumber );
	} else if (this->contains_scalar(key)) {
		return this->make_scalar_header_line( key );
	} else {
		std::ostringstream oss;
		oss << "No property " << key << " in atmosphere!";
		throw std::invalid_argument( oss.str() );
	}
}

std::string NCPA::Atmosphere1D::format_header_line( size_t column,
		const std::string &key, const std::string &units ) {
	std::ostringstream oss;
	oss 	<< "#% " << column
			<< ", " << key
			<< ", " << units;
	return oss.str();
}

std::string NCPA::Atmosphere1D::format_header_line( size_t column,
		const std::string &key, const std::string &units,
		double val ) {
	std::ostringstream oss;
	oss << NCPA::Atmosphere1D::format_header_line( column, key, units )
		<< ", " << NCPA::Atmosphere1D::format_for_stream(val);
	return oss.str();
}

std::string NCPA::Atmosphere1D::format_header_line( size_t column,
		const std::string &key, const std::string &units,
		int val ) {
	std::ostringstream oss;
	oss << NCPA::Atmosphere1D::format_header_line( column, key, units )
		<< ", " << NCPA::Atmosphere1D::format_for_stream(val);
	return oss.str();
}

std::string NCPA::Atmosphere1D::format_for_stream( double val ) {
	std::ostringstream os;
	if (std::fabs(std::log10(val)) >= NCPA_ATMOSPHERE_FIXED_MAX_LOG) {
		os 		<< std::scientific
				<< std::setprecision(NCPA_ATMOSPHERE_SCIENTIFIC_PRECISION);
	} else {
		os 	<< std::fixed
			<< std::setprecision(NCPA_ATMOSPHERE_FIXED_PRECISION);
	}
	os << val;
	return os.str();
}

std::string NCPA::Atmosphere1D::format_for_stream( int val ) {
	std::ostringstream os;
	if (std::fabs(std::log10(val)) >= NCPA_ATMOSPHERE_FIXED_MAX_LOG) {
		os 		<< std::scientific
				<< std::setprecision(NCPA_ATMOSPHERE_SCIENTIFIC_PRECISION)
				<< (double)val;
	} else {
		os 	<< std::fixed
			<< std::setprecision(NCPA_ATMOSPHERE_FIXED_PRECISION)
			<< val;
	}
	return os.str();
}
