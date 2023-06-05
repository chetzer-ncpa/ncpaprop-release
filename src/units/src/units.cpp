#include "units.h"
#include "util.h"
#include <map>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <iomanip>

// initialize static members
conversion_map_t NCPA::Units::map_ = conversion_map_t();
units_to_string_map_t NCPA::Units::enum_to_string_map_ = units_to_string_map_t();
units_to_string_map_t NCPA::Units::enum_to_abbr_map_ = units_to_string_map_t();
string_to_units_map_t NCPA::Units::string_to_enum_map_ = string_to_units_map_t();


void NCPA::Units::list_recognized_strings( std::ostream& os ) {
	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	int maxwidth = 0;
	string_to_units_map_t::const_iterator i;
	for ( i = NCPA::Units::string_to_enum_map_.cbegin(); i != NCPA::Units::string_to_enum_map_.cend(); ++i ) {
		maxwidth = NCPA::max<int>( maxwidth, (int)(*i).first.size() );
	}

	os << "Note: strings are not case-sensitive" << std::endl;
	os << std::setw( maxwidth ) << std::right << "String" << " : " << "Units" << std::endl;
	os << std::setw( maxwidth ) << std::right << "------" << " : " << "-----" << std::endl;
	for ( i = NCPA::Units::string_to_enum_map_.cbegin(); i != NCPA::Units::string_to_enum_map_.cend(); ++i ) {
		os << std::setw( maxwidth ) << std::right << (*i).first << " | " << NCPA::Units::toString( (*i).second ) << std::endl;
	}
}

std::vector<std::string> NCPA::Units::get_recognized_strings() {
	std::vector<std::string> names( NCPA::Units::string_to_enum_map_.size() );
	string_to_units_map_t::const_iterator i;
	for ( i = NCPA::Units::string_to_enum_map_.cbegin();
			i != NCPA::Units::string_to_enum_map_.cend(); ++i ) {
		names.push_back( i->first );
	}
	return names;
}

NCPA::units_t NCPA::Units::fromString( std::string s ) {

	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	if (NCPA::Units::string_to_enum_map_.count( NCPA::toUpperCase( s ) ) == 1) {
		return NCPA::Units::string_to_enum_map_.at( NCPA::toUpperCase( s ) );
	} else {
		return NCPA::UNITS_NONE;
		//throw std::out_of_range( "Unrecognized units string: " + s );
	}
	/*
	try {
		return NCPA::Units::string_to_enum_map_.at( NCPA::toUpperCase( s ) );
	} catch (std::out_of_range e) {
		throw std::out_of_range( "Unrecognized units string: " + s );
	}
	*/
}

std::string NCPA::Units::toString( units_t u ) {

	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	if (NCPA::Units::enum_to_string_map_.count( u ) == 1) {
		return NCPA::Units::enum_to_string_map_.at( u );
	} else {
		return "";
		//throw std::out_of_range( "Unrecognized units string: " + s );
	}
/*
	try {
		return NCPA::Units::enum_to_string_map_.at( u );
	} catch (std::out_of_range e) {
		throw std::out_of_range( "Unrecognized units type" );
	}
*/

}


std::string NCPA::Units::toStr( units_t u ) {

	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	if (NCPA::Units::enum_to_abbr_map_.count( u ) == 1) {
		return NCPA::Units::enum_to_abbr_map_.at( u );
	} else {
		return "";
		//throw std::out_of_range( "Unrecognized units string: " + s );
	}
/*
	try {
		return NCPA::Units::enum_to_abbr_map_.at( u );
	} catch (std::out_of_range e) {
		throw std::out_of_range( "Unrecognized units type" );
	}
*/
}


/*
Constructor for the UnitConverter class.
Internally the class uses a std::map to associate a
std::pair< units_t, units_t > key with a function pointer that performs the
indicated conversion.  The constructor populates this map with all defined
conversions.  No memory is dynamically allocated within the constructor so no
explicit destructor is required.  Also initializes string <---> enum maps.
*/
void NCPA::Units::initialize_() {

	map_.clear();
	// map_d1_.clear();
	// map_d2_.clear();
	enum_to_string_map_.clear();
	enum_to_abbr_map_.clear();
	string_to_enum_map_.clear();

	// conversion functions.  Functions must take in one double, put out one double, and not depend on any external variables.
	map_ = {
		{ get_unit_pair_( UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_FAHRENHEIT ), []( double in ) { return ( in * 1.8 ) + 32.0; } },
		{ get_unit_pair_( UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_KELVIN ), []( double in ) { return in + 273.15; } },
		{ get_unit_pair_( UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_FAHRENHEIT ), []( double in ) { return ( ( in - 273.15 ) * 1.8 ) + 32.0; } },
		{ get_unit_pair_( UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_CELSIUS ), []( double in ) { return in - 273.15; } },
		{ get_unit_pair_( UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_KELVIN ), []( double in ) { return ( ( in - 32.0 ) * 5.0 / 9.0 ) + 273.15; } },
		{ get_unit_pair_( UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_CELSIUS ), []( double in ) { return ( in - 32.0 ) * 5.0 / 9.0; } },
		{ get_unit_pair_( UNITS_DISTANCE_METERS, UNITS_DISTANCE_KILOMETERS ), []( double in ) { return in * 0.001; } },
		{ get_unit_pair_( UNITS_DISTANCE_KILOMETERS, UNITS_DISTANCE_METERS ), []( double in ) { return in * 1000.0; } },
		{ get_unit_pair_( UNITS_SPEED_METERS_PER_SECOND, UNITS_SPEED_KILOMETERS_PER_SECOND ), []( double in ) { return in * 0.001; } },
		{ get_unit_pair_( UNITS_SPEED_KILOMETERS_PER_SECOND, UNITS_SPEED_METERS_PER_SECOND ), []( double in ) { return in * 1000.0; } },
		{ get_unit_pair_( UNITS_PRESSURE_PASCALS, UNITS_PRESSURE_MILLIBARS ), []( double in ) { return in * 0.01; } },
		{ get_unit_pair_( UNITS_PRESSURE_PASCALS, UNITS_PRESSURE_HECTOPASCALS ), []( double in ) { return in * 0.01; } },
		{ get_unit_pair_( UNITS_PRESSURE_MILLIBARS, UNITS_PRESSURE_PASCALS ), []( double in ) { return in * 100.0; } },
		{ get_unit_pair_( UNITS_PRESSURE_MILLIBARS, UNITS_PRESSURE_HECTOPASCALS ), []( double in ) { return in; } },
		{ get_unit_pair_( UNITS_PRESSURE_HECTOPASCALS, UNITS_PRESSURE_PASCALS ), []( double in ) { return in * 100; } },
		{ get_unit_pair_( UNITS_PRESSURE_HECTOPASCALS, UNITS_PRESSURE_MILLIBARS ), []( double in ) { return in; } },
		{ get_unit_pair_( UNITS_PRESSURE_ATMOSPHERES, UNITS_PRESSURE_HECTOPASCALS ), []( double in ) { return in * 1013.25; } },
		{ get_unit_pair_( UNITS_PRESSURE_ATMOSPHERES, UNITS_PRESSURE_PASCALS ), []( double in ) { return in * 101325.0; } },
		{ get_unit_pair_( UNITS_PRESSURE_ATMOSPHERES, UNITS_PRESSURE_MILLIBARS ), []( double in ) { return in * 1013.25; } },
		{ get_unit_pair_( UNITS_PRESSURE_HECTOPASCALS, UNITS_PRESSURE_ATMOSPHERES ), []( double in ) { return in * 0.000986923; } },
		{ get_unit_pair_( UNITS_PRESSURE_PASCALS, UNITS_PRESSURE_ATMOSPHERES ), []( double in ) { return in * 0.00000986923; } },
		{ get_unit_pair_( UNITS_PRESSURE_MILLIBARS, UNITS_PRESSURE_ATMOSPHERES ), []( double in ) { return in * 0.000986923; } },
		{ get_unit_pair_( UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER, UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER ), []( double in ) { return in * 0.001; } },
		{ get_unit_pair_( UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER, UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER ), []( double in ) { return in * 1000.0; } },
		{ get_unit_pair_( UNITS_ANGLE_RADIANS, UNITS_ANGLE_DEGREES ), []( double in ) { return in * 180.0 / PI; } },
		{ get_unit_pair_( UNITS_ANGLE_DEGREES, UNITS_ANGLE_RADIANS ), []( double in ) { return in * PI / 180.0; } },
		{ get_unit_pair_( UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST ), []( double in ) {
			double out = 90.0 - in;
			while (out < 0) {
				out += 360.0;
			}
			while (out >= 360.0) {
				out -= 360.0;
			}
			return out;
		} },
		{ get_unit_pair_( UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH ), []( double in ) {
			double out = 90.0 - in;
			while (out < 0) {
				out += 360.0;
			}
			while (out >= 360.0) {
				out -= 360.0;
			}
			return out;
		} }
	};



	// take an enum value, output it's string value (full version)
	enum_to_string_map_ = {
		{ UNITS_NONE, "" },
		{ UNITS_TEMPERATURE_KELVIN, "degrees Kelvin" },
		{ UNITS_TEMPERATURE_CELSIUS, "degrees Celsius" },
		{ UNITS_TEMPERATURE_FAHRENHEIT, "degrees Fahrenheit" },
		{ UNITS_DISTANCE_METERS, "meters" },
		{ UNITS_DISTANCE_KILOMETERS, "kilometers" },
		{ UNITS_SPEED_METERS_PER_SECOND, "meters per second" },
		{ UNITS_SPEED_KILOMETERS_PER_SECOND, "kilometers per second" },
		{ UNITS_PRESSURE_PASCALS, "Pascals" },
		{ UNITS_PRESSURE_MILLIBARS, "millibars" },
		{ UNITS_PRESSURE_HECTOPASCALS, "hectopascals" },
		{ UNITS_PRESSURE_ATMOSPHERES, "atmospheres" },
		{ UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER, "kilograms per cubic meter" },
		{ UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER, "grams per cubic centimeter" },
		{ UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, "degrees clockwise from North" },
		{ UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, "degrees counterclockwise from East" },
		{ UNITS_ANGLE_DEGREES, "degrees" },
		{ UNITS_ANGLE_RADIANS, "radians" },
	};

	enum_to_abbr_map_ = {
		{ UNITS_NONE, "" },
		{ UNITS_TEMPERATURE_KELVIN, "K" },
		{ UNITS_TEMPERATURE_CELSIUS, "C" },
		{ UNITS_TEMPERATURE_FAHRENHEIT, "F" },
		{ UNITS_DISTANCE_METERS, "m" },
		{ UNITS_DISTANCE_KILOMETERS, "km" },
		{ UNITS_SPEED_METERS_PER_SECOND, "m/s" },
		{ UNITS_SPEED_KILOMETERS_PER_SECOND, "km/s" },
		{ UNITS_PRESSURE_PASCALS, "Pa" },
		{ UNITS_PRESSURE_MILLIBARS, "mbar" },
		{ UNITS_PRESSURE_HECTOPASCALS, "hPa" },
		{ UNITS_PRESSURE_ATMOSPHERES, "atm" },
		{ UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER, "kg/m3" },
		{ UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER, "g/cm3" },
		{ UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, "deg CW from N" },
		{ UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, "deg CCW from E" },
		{ UNITS_ANGLE_DEGREES, "deg" },
		{ UNITS_ANGLE_RADIANS, "rad" },
	};

	string_to_enum_map_ = {
		{ "", UNITS_NONE },
		{ "N/A", UNITS_NONE },
		{ "K", UNITS_TEMPERATURE_KELVIN },
		{ "DEGK", UNITS_TEMPERATURE_KELVIN },
		{ "DEG K", UNITS_TEMPERATURE_KELVIN },
		{ "DEGREES K", UNITS_TEMPERATURE_KELVIN },
		{ "KELVIN", UNITS_TEMPERATURE_KELVIN },
		{ "C", UNITS_TEMPERATURE_CELSIUS },
		{ "DEGC", UNITS_TEMPERATURE_CELSIUS },
		{ "DEG C", UNITS_TEMPERATURE_CELSIUS },
		{ "DEGREES C", UNITS_TEMPERATURE_CELSIUS },
		{ "CELSIUS", UNITS_TEMPERATURE_CELSIUS },
		{ "CENTIGRADE", UNITS_TEMPERATURE_CELSIUS },
		{ "F", UNITS_TEMPERATURE_FAHRENHEIT },
		{ "DEGF", UNITS_TEMPERATURE_FAHRENHEIT },
		{ "DEG F", UNITS_TEMPERATURE_FAHRENHEIT },
		{ "DEGREES F", UNITS_TEMPERATURE_FAHRENHEIT },
		{ "FAHRENHEIT", UNITS_TEMPERATURE_FAHRENHEIT },
		{ "M", UNITS_DISTANCE_METERS },
		{ "METERS", UNITS_DISTANCE_METERS },
		{ "KM", UNITS_DISTANCE_KILOMETERS },
		{ "KILOMETERS", UNITS_DISTANCE_KILOMETERS },
		{ "M/S", UNITS_SPEED_METERS_PER_SECOND },
		{ "MPS", UNITS_SPEED_METERS_PER_SECOND },
		{ "MPERS", UNITS_SPEED_METERS_PER_SECOND },
		{ "M PER S", UNITS_SPEED_METERS_PER_SECOND },
		{ "METERS PER SECOND", UNITS_SPEED_METERS_PER_SECOND },
		{ "KM/S", UNITS_SPEED_KILOMETERS_PER_SECOND },
		{ "KMPS", UNITS_SPEED_KILOMETERS_PER_SECOND },
		{ "KMPERS", UNITS_SPEED_KILOMETERS_PER_SECOND },
		{ "KM PER S", UNITS_SPEED_KILOMETERS_PER_SECOND },
		{ "KILOMETERS PER SECOND", UNITS_SPEED_KILOMETERS_PER_SECOND },
		{ "PA", UNITS_PRESSURE_PASCALS },
		{ "PASCAL", UNITS_PRESSURE_PASCALS },
		{ "PASCALS", UNITS_PRESSURE_PASCALS },
		{ "MBAR", UNITS_PRESSURE_MILLIBARS },
		{ "MILLIBAR", UNITS_PRESSURE_MILLIBARS },
		{ "MILLIBARS", UNITS_PRESSURE_MILLIBARS },
		{ "HECTOPASCALS", UNITS_PRESSURE_HECTOPASCALS },
		{ "HPA", UNITS_PRESSURE_HECTOPASCALS },
		{ "ATMOSPHERES", UNITS_PRESSURE_ATMOSPHERES },
		{ "ATM", UNITS_PRESSURE_ATMOSPHERES },
		{ "KG/M3", UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER },
		{ "KGPM3", UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER },
		{ "KILOGRAMS PER CUBIC METER", UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER },
		{ "G/CM3", UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER },
		{ "GPCM3", UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER },
		{ "GRAMS PER CUBIC CENTIMETER", UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER },
		{ "DEGREES CLOCKWISE FROM NORTH", UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH },
		{ "DEG CW FROM N", UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH },
		{ "AZIMUTH", UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH },
		{ "DEGREES COUNTERCLOCKWISE FROM EAST", UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST },
		{ "DEG CCW FROM E", UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST },
		{ "DEG", UNITS_ANGLE_DEGREES },
		{ "DEGREES", UNITS_ANGLE_DEGREES },
		{ "RAD", UNITS_ANGLE_RADIANS },
		{ "RADIANS", UNITS_ANGLE_RADIANS }
	};

}

bool NCPA::Units::ready_() {
	return !(map_.empty());
}

double NCPA::Units::convert( double in, NCPA::units_t type_in, NCPA::units_t type_out ) {
	double out = 0.0;
	NCPA::Units::convert( &in, 1, type_in, type_out, &out );
	/*
	try {
		// Call the vector conversion with one-element vectors
		NCPA::Units::convert( &in, 1, type_in, type_out, &out );
	} catch (const std::out_of_range& oor) {
		// didn't find the requested conversion, kick it upstairs, user
		// will have been notified in the called function
		throw;
	}
	*/
	return out;
}

/*
Converts one or more double values from one unit to another.
*/
void NCPA::Units::convert( const double *in, unsigned int nSamples,
	NCPA::units_t type_in, NCPA::units_t type_out, double *out ) {

	// Is the conversion from one type to the same type?
	if (type_in == type_out) {

		// If the identity conversion is in place, just return
		if (in == out) {
			return;
		}

		// Copy the old values to the new values with no conversion
		std::memcpy( out, in, nSamples*sizeof(double) );
		return;
	}

	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	// Create a pair with the requested in and out units
	conversion_pair cpair( type_in, type_out );

	// function will go here if found
	conversion_function fPtr;
	try {
		// see if the requested conversion exists in the map.  If not, it throws
		// an out_of_range exception that is caught later
		fPtr = NCPA::Units::map_.at( cpair );

		// perform the requested conversion
		for (unsigned int i = 0; i < nSamples; i++) {
			out[ i ] = fPtr( in[ i ] );
		}

	} catch (const std::out_of_range& oor) {
		// didn't find the conversion, notify the user and kick it upstairs
		throw std::out_of_range( "Undefined conversion requested from "
			+ NCPA::Units::toString( type_in ) + " to " + NCPA::Units::toString( type_out ) );
	}
}


/*
Generates and returns a std::pair with the two unit types, for use as a map key.
*/
conversion_pair NCPA::Units::get_unit_pair_( NCPA::units_t t1, NCPA::units_t t2 ) {
	return std::make_pair( t1, t2 );
}






