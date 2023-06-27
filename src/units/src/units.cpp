#include "NCPACommon.h"
#include "units.h"
#include <map>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>

// initialize static members
conversion_map_t NCPA::Units::map_ = conversion_map_t();
units_to_string_map_t NCPA::Units::enum_to_string_map_ = units_to_string_map_t();
units_to_string_map_t NCPA::Units::enum_to_abbr_map_ = units_to_string_map_t();
string_to_units_map_t NCPA::Units::string_to_enum_map_ = string_to_units_map_t();

NCPA::invalid_conversion::invalid_conversion( NCPA::units_t from, NCPA::units_t to )
	: std::out_of_range("Invalid conversion from " + NCPA::Units::toString(from)
						+ " to " + NCPA::Units::toString(to)) {}

NCPA::invalid_conversion::invalid_conversion( const std::string &msg )
	: std::out_of_range(msg) {}

NCPA::invalid_conversion::invalid_conversion() : std::out_of_range("Invalid conversion") {}

bool operator==(NCPA::units_t a, NCPA::units_t b) {
	return (static_cast<size_t>(a) == static_cast<size_t>(b));
}

bool operator!=(NCPA::units_t a, NCPA::units_t b) {
	return !(a == b);
}

bool NCPA::equal_or_none( NCPA::units_t a, NCPA::units_t b ) {
	return (a == b) || (a == NCPA::units_t::NONE) || (b == NCPA::units_t::NONE);
}

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
		return NCPA::units_t::NONE;
	}
}

std::string NCPA::Units::toString( units_t u ) {

	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	if (NCPA::Units::enum_to_string_map_.count( u ) == 1) {
		return NCPA::Units::enum_to_string_map_.at( u );
	} else {
		return "";
	}
}


std::string NCPA::Units::toStr( units_t u ) {

	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}

	if (NCPA::Units::enum_to_abbr_map_.count( u ) == 1) {
		return NCPA::Units::enum_to_abbr_map_.at( u );
	} else {
		return "";
	}
}

NCPA::units_t NCPA::Units::null() {
	return NCPA::units_t::NONE;
}


/*
Constructor for the Units class.
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
		{ get_unit_pair_( NCPA::units_t::TEMPERATURE_CELSIUS, NCPA::units_t::TEMPERATURE_FAHRENHEIT ), []( double in ) { return ( in * 1.8 ) + 32.0; } },
		{ get_unit_pair_( NCPA::units_t::TEMPERATURE_CELSIUS, NCPA::units_t::TEMPERATURE_KELVIN ), []( double in ) { return in + 273.15; } },
		{ get_unit_pair_( NCPA::units_t::TEMPERATURE_KELVIN, NCPA::units_t::TEMPERATURE_FAHRENHEIT ), []( double in ) { return ( ( in - 273.15 ) * 1.8 ) + 32.0; } },
		{ get_unit_pair_( NCPA::units_t::TEMPERATURE_KELVIN, NCPA::units_t::TEMPERATURE_CELSIUS ), []( double in ) { return in - 273.15; } },
		{ get_unit_pair_( NCPA::units_t::TEMPERATURE_FAHRENHEIT, NCPA::units_t::TEMPERATURE_KELVIN ), []( double in ) { return ( ( in - 32.0 ) * 5.0 / 9.0 ) + 273.15; } },
		{ get_unit_pair_( NCPA::units_t::TEMPERATURE_FAHRENHEIT, NCPA::units_t::TEMPERATURE_CELSIUS ), []( double in ) { return ( in - 32.0 ) * 5.0 / 9.0; } },
		{ get_unit_pair_( NCPA::units_t::DISTANCE_METERS, NCPA::units_t::DISTANCE_KILOMETERS ), []( double in ) { return in * 0.001; } },
		{ get_unit_pair_( NCPA::units_t::DISTANCE_KILOMETERS, NCPA::units_t::DISTANCE_METERS ), []( double in ) { return in * 1000.0; } },
		{ get_unit_pair_( NCPA::units_t::SPEED_METERS_PER_SECOND, NCPA::units_t::SPEED_KILOMETERS_PER_SECOND ), []( double in ) { return in * 0.001; } },
		{ get_unit_pair_( NCPA::units_t::SPEED_KILOMETERS_PER_SECOND, NCPA::units_t::SPEED_METERS_PER_SECOND ), []( double in ) { return in * 1000.0; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_PASCALS, NCPA::units_t::PRESSURE_MILLIBARS ), []( double in ) { return in * 0.01; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_PASCALS, NCPA::units_t::PRESSURE_HECTOPASCALS ), []( double in ) { return in * 0.01; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_MILLIBARS, NCPA::units_t::PRESSURE_PASCALS ), []( double in ) { return in * 100.0; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_MILLIBARS, NCPA::units_t::PRESSURE_HECTOPASCALS ), []( double in ) { return in; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_HECTOPASCALS, NCPA::units_t::PRESSURE_PASCALS ), []( double in ) { return in * 100; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_HECTOPASCALS, NCPA::units_t::PRESSURE_MILLIBARS ), []( double in ) { return in; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_ATMOSPHERES, NCPA::units_t::PRESSURE_HECTOPASCALS ), []( double in ) { return in * 1013.25; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_ATMOSPHERES, NCPA::units_t::PRESSURE_PASCALS ), []( double in ) { return in * 101325.0; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_ATMOSPHERES, NCPA::units_t::PRESSURE_MILLIBARS ), []( double in ) { return in * 1013.25; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_HECTOPASCALS, NCPA::units_t::PRESSURE_ATMOSPHERES ), []( double in ) { return in * 0.000986923; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_PASCALS, NCPA::units_t::PRESSURE_ATMOSPHERES ), []( double in ) { return in * 0.00000986923; } },
		{ get_unit_pair_( NCPA::units_t::PRESSURE_MILLIBARS, NCPA::units_t::PRESSURE_ATMOSPHERES ), []( double in ) { return in * 0.000986923; } },
		{ get_unit_pair_( NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER, NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER ), []( double in ) { return in * 0.001; } },
		{ get_unit_pair_( NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER, NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER ), []( double in ) { return in * 1000.0; } },
		{ get_unit_pair_( NCPA::units_t::ANGLE_RADIANS, NCPA::units_t::ANGLE_DEGREES ), []( double in ) { return in * 180.0 / M_PI; } },
		{ get_unit_pair_( NCPA::units_t::ANGLE_DEGREES, NCPA::units_t::ANGLE_RADIANS ), []( double in ) { return in * M_PI / 180.0; } },
		{ get_unit_pair_( NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST ), []( double in ) {
			double out = 90.0 - in;
			while (out < 0) {
				out += 360.0;
			}
			while (out >= 360.0) {
				out -= 360.0;
			}
			return out;
		} },
		{ get_unit_pair_( NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH ), []( double in ) {
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
		{ NCPA::units_t::NONE, "" },
		{ NCPA::units_t::TEMPERATURE_KELVIN, "degrees Kelvin" },
		{ NCPA::units_t::TEMPERATURE_CELSIUS, "degrees Celsius" },
		{ NCPA::units_t::TEMPERATURE_FAHRENHEIT, "degrees Fahrenheit" },
		{ NCPA::units_t::DISTANCE_METERS, "meters" },
		{ NCPA::units_t::DISTANCE_KILOMETERS, "kilometers" },
		{ NCPA::units_t::SPEED_METERS_PER_SECOND, "meters per second" },
		{ NCPA::units_t::SPEED_KILOMETERS_PER_SECOND, "kilometers per second" },
		{ NCPA::units_t::PRESSURE_PASCALS, "Pascals" },
		{ NCPA::units_t::PRESSURE_MILLIBARS, "millibars" },
		{ NCPA::units_t::PRESSURE_HECTOPASCALS, "hectopascals" },
		{ NCPA::units_t::PRESSURE_ATMOSPHERES, "atmospheres" },
		{ NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER, "kilograms per cubic meter" },
		{ NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER, "grams per cubic centimeter" },
		{ NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, "degrees clockwise from North" },
		{ NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, "degrees counterclockwise from East" },
		{ NCPA::units_t::ANGLE_DEGREES, "degrees" },
		{ NCPA::units_t::ANGLE_RADIANS, "radians" },
	};

	enum_to_abbr_map_ = {
		{ NCPA::units_t::NONE, "" },
		{ NCPA::units_t::TEMPERATURE_KELVIN, "K" },
		{ NCPA::units_t::TEMPERATURE_CELSIUS, "C" },
		{ NCPA::units_t::TEMPERATURE_FAHRENHEIT, "F" },
		{ NCPA::units_t::DISTANCE_METERS, "m" },
		{ NCPA::units_t::DISTANCE_KILOMETERS, "km" },
		{ NCPA::units_t::SPEED_METERS_PER_SECOND, "m/s" },
		{ NCPA::units_t::SPEED_KILOMETERS_PER_SECOND, "km/s" },
		{ NCPA::units_t::PRESSURE_PASCALS, "Pa" },
		{ NCPA::units_t::PRESSURE_MILLIBARS, "mbar" },
		{ NCPA::units_t::PRESSURE_HECTOPASCALS, "hPa" },
		{ NCPA::units_t::PRESSURE_ATMOSPHERES, "atm" },
		{ NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER, "kg/m3" },
		{ NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER, "g/cm3" },
		{ NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, "deg CW from N" },
		{ NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, "deg CCW from E" },
		{ NCPA::units_t::ANGLE_DEGREES, "deg" },
		{ NCPA::units_t::ANGLE_RADIANS, "rad" },
	};

	string_to_enum_map_ = {
		{ "", NCPA::units_t::NONE },
		{ "N/A", NCPA::units_t::NONE },
		{ "K", NCPA::units_t::TEMPERATURE_KELVIN },
		{ "DEGK", NCPA::units_t::TEMPERATURE_KELVIN },
		{ "DEG K", NCPA::units_t::TEMPERATURE_KELVIN },
		{ "DEGREES K", NCPA::units_t::TEMPERATURE_KELVIN },
		{ "KELVIN", NCPA::units_t::TEMPERATURE_KELVIN },
		{ "C", NCPA::units_t::TEMPERATURE_CELSIUS },
		{ "DEGC", NCPA::units_t::TEMPERATURE_CELSIUS },
		{ "DEG C", NCPA::units_t::TEMPERATURE_CELSIUS },
		{ "DEGREES C", NCPA::units_t::TEMPERATURE_CELSIUS },
		{ "CELSIUS", NCPA::units_t::TEMPERATURE_CELSIUS },
		{ "CENTIGRADE", NCPA::units_t::TEMPERATURE_CELSIUS },
		{ "F", NCPA::units_t::TEMPERATURE_FAHRENHEIT },
		{ "DEGF", NCPA::units_t::TEMPERATURE_FAHRENHEIT },
		{ "DEG F", NCPA::units_t::TEMPERATURE_FAHRENHEIT },
		{ "DEGREES F", NCPA::units_t::TEMPERATURE_FAHRENHEIT },
		{ "FAHRENHEIT", NCPA::units_t::TEMPERATURE_FAHRENHEIT },
		{ "M", NCPA::units_t::DISTANCE_METERS },
		{ "METERS", NCPA::units_t::DISTANCE_METERS },
		{ "KM", NCPA::units_t::DISTANCE_KILOMETERS },
		{ "KILOMETERS", NCPA::units_t::DISTANCE_KILOMETERS },
		{ "M/S", NCPA::units_t::SPEED_METERS_PER_SECOND },
		{ "MPS", NCPA::units_t::SPEED_METERS_PER_SECOND },
		{ "MPERS", NCPA::units_t::SPEED_METERS_PER_SECOND },
		{ "M PER S", NCPA::units_t::SPEED_METERS_PER_SECOND },
		{ "METERS PER SECOND", NCPA::units_t::SPEED_METERS_PER_SECOND },
		{ "KM/S", NCPA::units_t::SPEED_KILOMETERS_PER_SECOND },
		{ "KMPS", NCPA::units_t::SPEED_KILOMETERS_PER_SECOND },
		{ "KMPERS", NCPA::units_t::SPEED_KILOMETERS_PER_SECOND },
		{ "KM PER S", NCPA::units_t::SPEED_KILOMETERS_PER_SECOND },
		{ "KILOMETERS PER SECOND", NCPA::units_t::SPEED_KILOMETERS_PER_SECOND },
		{ "PA", NCPA::units_t::PRESSURE_PASCALS },
		{ "PASCAL", NCPA::units_t::PRESSURE_PASCALS },
		{ "PASCALS", NCPA::units_t::PRESSURE_PASCALS },
		{ "MBAR", NCPA::units_t::PRESSURE_MILLIBARS },
		{ "MILLIBAR", NCPA::units_t::PRESSURE_MILLIBARS },
		{ "MILLIBARS", NCPA::units_t::PRESSURE_MILLIBARS },
		{ "HECTOPASCALS", NCPA::units_t::PRESSURE_HECTOPASCALS },
		{ "HPA", NCPA::units_t::PRESSURE_HECTOPASCALS },
		{ "ATMOSPHERES", NCPA::units_t::PRESSURE_ATMOSPHERES },
		{ "ATM", NCPA::units_t::PRESSURE_ATMOSPHERES },
		{ "KG/M3", NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER },
		{ "KGPM3", NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER },
		{ "KILOGRAMS PER CUBIC METER", NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER },
		{ "G/CM3", NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER },
		{ "GPCM3", NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER },
		{ "GRAMS PER CUBIC CENTIMETER", NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER },
		{ "DEGREES CLOCKWISE FROM NORTH", NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH },
		{ "DEG CW FROM N", NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH },
		{ "AZIMUTH", NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH },
		{ "DEGREES COUNTERCLOCKWISE FROM EAST", NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST },
		{ "DEG CCW FROM E", NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST },
		{ "DEG", NCPA::units_t::ANGLE_DEGREES },
		{ "DEGREES", NCPA::units_t::ANGLE_DEGREES },
		{ "RAD", NCPA::units_t::ANGLE_RADIANS },
		{ "RADIANS", NCPA::units_t::ANGLE_RADIANS }
	};

}

bool NCPA::Units::ready_() {
	return !(map_.empty());
}

double NCPA::Units::convert( double in, NCPA::units_t type_in, const std::string &type_out ) {
	return NCPA::Units::convert(
			in, type_in, NCPA::Units::fromString(type_out) );
}

double NCPA::Units::convert( double in, const std::string &type_in, NCPA::units_t type_out ) {
	return NCPA::Units::convert(
			in, NCPA::Units::fromString(type_in), type_out );
}

double NCPA::Units::convert( double in, const std::string &type_in, const std::string &type_out ) {
	return NCPA::Units::convert(
			in, NCPA::Units::fromString(type_in), NCPA::Units::fromString(type_out) );
}

double NCPA::Units::convert( double in, NCPA::units_t type_in, NCPA::units_t type_out ) {
	double out = 0.0;
	NCPA::Units::convert( &in, 1, type_in, type_out, &out );
	return out;
}

void NCPA::Units::convert( const double *in, unsigned int nSamples,
	const std::string &type_in, const std::string &type_out, double *out ) {
	NCPA::Units::convert(
			in, nSamples, NCPA::Units::fromString(type_in),
			NCPA::Units::fromString(type_out), out
			);
}

void NCPA::Units::convert( const double *in, unsigned int nSamples,
	NCPA::units_t type_in, const std::string &type_out, double *out ) {
	NCPA::Units::convert(
			in, nSamples, type_in,
			NCPA::Units::fromString(type_out), out
			);
}

void NCPA::Units::convert( const double *in, unsigned int nSamples,
	const std::string &type_in, NCPA::units_t type_out, double *out ) {
	NCPA::Units::convert(
			in, nSamples, NCPA::Units::fromString(type_in),
			type_out, out
			);
}

/*
Converts one or more double values from one unit to another.
*/
void NCPA::Units::convert( const double *in, unsigned int nSamples,
	NCPA::units_t type_in, NCPA::units_t type_out, double *out ) {

	// Is the conversion from one type to the same type, or to/from NONE?
	if (NCPA::equal_or_none(type_in,type_out)) {

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
		throw NCPA::invalid_conversion(type_in,type_out);
	}
}

bool NCPA::Units::can_convert( NCPA::units_t type_in, NCPA::units_t type_out ) {
	// Is the conversion from one type to the same type?
	if (NCPA::equal_or_none(type_in,type_out)) {
		return true;
	}
	conversion_pair cpair( type_in, type_out );
	return (NCPA::Units::map_.find(cpair) != NCPA::Units::map_.end());
}

bool NCPA::Units::can_convert( const std::string &type_in, NCPA::units_t type_out ) {
	return NCPA::Units::can_convert(
			NCPA::Units::fromString(type_in),
			type_out );
}

bool NCPA::Units::can_convert( NCPA::units_t type_in, const std::string &type_out ) {
	return NCPA::Units::can_convert(
			type_in,
			NCPA::Units::fromString(type_out) );
}

bool NCPA::Units::can_convert( const std::string &type_in, const std::string &type_out ) {
	return NCPA::Units::can_convert(
			NCPA::Units::fromString(type_in),
			NCPA::Units::fromString(type_out) );
}


/*
Generates and returns a std::pair with the two unit types, for use as a map key.
*/
conversion_pair NCPA::Units::get_unit_pair_( NCPA::units_t t1, NCPA::units_t t2 ) {
	return std::make_pair( t1, t2 );
}






