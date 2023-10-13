/*
Classes, functions and constants for unit conversion.

Classes:
NCPA::UnitConverter: 	Handle all defined unit conversions.  Can convert single values or vectors.
			Will throw an exception if an undefined conversion is requested.

Constants:
NCPA::units_t	NCPA::units_t::NONE
NCPA::units_t 	NCPA::units_t::TEMPERATURE_CELSIUS
NCPA::units_t 	NCPA::units_t::TEMPERATURE_KELVIN
NCPA::units_t 	NCPA::units_t::TEMPERATURE_FAHRENHEIT
NCPA::units_t 	NCPA::units_t::DISTANCE_METERS
NCPA::units_t 	NCPA::units_t::DISTANCE_KILOMETERS
NCPA::units_t 	NCPA::units_t::SPEED_METERS_PER_SECOND
NCPA::units_t 	NCPA::units_t::SPEED_KILOMETERS_PER_SECOND
NCPA::units_t 	NCPA::units_t::PRESSURE_PASCALS
NCPA::units_t 	NCPA::units_t::PRESSURE_MILLIBARS
NCPA::units_t   NCPA::units_t::PRESSURE_HECTOPASCALS
NCPA::units_t 	NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER
NCPA::units_t 	NCPA::units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER
NCPA::units_t	NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH
NCPA::units_t	NCPA::units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST
NCPA::units_t 	NCPA::units_t::ANGLE_DEGREES
NCPA::units_t 	NCPA::units_t::ANGLE_RADIANS
NCPA::units_t 	NCPA::units_t::ATTENUATION_NEPERS_PER_METER
NCPA::units_t 	NCPA::units_t::ATTENUATION_DECIBELS_PER_METER
NCPA::units_t 	NCPA::units_t::ATTENUATION_DECIBELS_PER_KILOMETER
NCPA::units_t	NCPA::units_t::RATIO_DECIMAL
NCPA::units_t	NCPA::units_t::RATIO_PERCENT

Examples:
	using namespace NCPA;

	// Convert a single value
	double temp_c = 50.0, temp_k;
	temp_k = Units::convert( temp_c, "C", "K" );
	// or
	Units::convert( &temp_c, 1, "C", "K", &temp_k );

	// Convert a vector of values
	double Pa_vec[ 20 ] = { ... };
	double mbar_vec[ 20 ];
	Units::convert( Pa_vec, 20, "Pa", "mbar", mbar_vec );

	// Can also convert in-place
	double distance[ 300 ] = { ... };
	Units::convert( distance, 300, "m", "km", distance );

	// if you need access to the actual type for some reason
	NCPA::units_t km = NCPA::units_t::DISTANCE_KILOMETERS;
	// or, easier to remember
	NCPA::units_t km = NCPA::Units::fromString("km");


To add a unit and its associated conversions, the following should be done:
1. Add symbol(s) to the units_t enum in this file
2. In units.cpp, add:
	* Unit conversion lambda functions to map_ variable
	* String to/from enum conversions to string_to_enum_map_, enum_to_string_map_, and
	  enum_to_abbr_map_ variables
*/

#ifndef NCPA_UNITS_H_INCLUDED
#define NCPA_UNITS_H_INCLUDED

#define HAVE_NCPA_UNITS_LIBRARY

#include <map>
#include <vector>
#include <stack>
#include <utility>
#include <iostream>
#include <stdexcept>

namespace NCPA {

	/**
	 * An enum of unit constants.
	 * Constants that can be used to identify or specify units.
	 */
	enum class units_t : size_t {
		NONE = 0,										/**< Indicates no units */

		TEMPERATURE_KELVIN,								/**< Temperature in Kelvin */
		TEMPERATURE_CELSIUS,							/**< Temperature in Celsius */
		TEMPERATURE_FAHRENHEIT,							/**< Temperature in Fahrenheit */

		DISTANCE_METERS,								/**< Distance in meters */
		DISTANCE_KILOMETERS,							/**< Distance in kilometers */

		SPEED_METERS_PER_SECOND,						/**< Speed in m/s */
		SPEED_KILOMETERS_PER_SECOND,					/**< Speed in km/s */

		PRESSURE_PASCALS,								/**< Pressure in Pa */
		PRESSURE_MILLIBARS,								/**< Pressure in mbar */
		PRESSURE_HECTOPASCALS,							/**< Pressure in hPa */
		PRESSURE_ATMOSPHERES,							/**< Pressure in atm */

		DENSITY_KILOGRAMS_PER_CUBIC_METER,				/**< Density in kg/m^3 */
		DENSITY_GRAMS_PER_CUBIC_CENTIMETER,				/**< Density in g/cm^3 */

		DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH,			/**< Direction in geographic azimuth */
		DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST,	/**< Direction in "math" convention */

		ANGLE_DEGREES,									/**< Angles in degrees */
		ANGLE_RADIANS,									/**< Angles in radians */

		ATTENUATION_NEPERS_PER_METER,					/**< Attenuation in np/m */
		ATTENUATION_DECIBELS_PER_METER,					/**< Attenuation in np/m */
		ATTENUATION_DECIBELS_PER_KILOMETER,				/**< Attenuation in dB/m */

		RATIO_DECIMAL,									/**< Ratio expressed as a decimal */
		RATIO_PERCENT,									/**< Ratio expressed as a percent */
	};

	class invalid_conversion : public std::out_of_range {
	public:
		invalid_conversion();
		invalid_conversion( const std::string &msg );
		invalid_conversion( units_t from, units_t to );
	};
}

#ifndef NULL_UNIT
#define NULL_UNIT NCPA::units_t::NONE
#endif

bool operator==(NCPA::units_t a, NCPA::units_t b);
bool operator!=(NCPA::units_t a, NCPA::units_t b);

typedef std::pair< NCPA::units_t, NCPA::units_t > conversion_pair;
typedef double (*conversion_function)(double);
typedef std::map< conversion_pair, conversion_function > conversion_map_t;
typedef std::map< std::string, NCPA::units_t > string_to_units_map_t;
typedef std::map< NCPA::units_t, std::string > units_to_string_map_t;

namespace NCPA {

	bool equal_or_none(units_t a, units_t b);

	/**
	 * A class for converting units.
	 */
	class Units {

	public:

		/**
		 * Returns whether a conversion for this pair has been defined.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return true if the conversion has been defined, false otherwise
		 */
		static bool can_convert( units_t type_in, units_t type_out );

		/**
		 * Returns whether a conversion for this pair has been defined.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return true if the conversion has been defined, false otherwise
		 */
		static bool can_convert( const std::string &type_in, units_t type_out );

		/**
		 * Returns whether a conversion for this pair has been defined.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return true if the conversion has been defined, false otherwise
		 */
		static bool can_convert( units_t type_in, const std::string &type_out );

		/**
		 * Returns whether a conversion for this pair has been defined.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return true if the conversion has been defined, false otherwise
		 */
		static bool can_convert( const std::string &type_in, const std::string &type_out );



		/**
		 * Convert an array of numbers from one unit to another.
		 * @param in 		A pointer to an array of double values
		 * @param nSamples 	The number of consecutive samples to convert
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @param out		A pointer to a preallocated array to hold the converted values
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 * @see units_t
		 */
		static void convert( const double *in, unsigned int nSamples,
			units_t type_in, units_t type_out, double *out );

		/**
		 * Convert an array of numbers from one unit to another.
		 * @param in 		A pointer to an array of double values
		 * @param nSamples 	The number of consecutive samples to convert
		 * @param type_in	String of the units to convert from
		 * @param type_out	String of the units to convert to
		 * @param out		A pointer to a preallocated array to hold the converted values
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 * @see units_t
		 */
		static void convert( const double *in, unsigned int nSamples,
			const std::string &type_in, const std::string &type_out, double *out );

		/**
		 * Convert an array of numbers from one unit to another.
		 * @param in 		A pointer to an array of double values
		 * @param nSamples 	The number of consecutive samples to convert
		 * @param type_in	The units to convert from
		 * @param type_out	String of the units to convert to
		 * @param out		A pointer to a preallocated array to hold the converted values
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 * @see units_t
		 */
		static void convert( const double *in, unsigned int nSamples,
			units_t type_in, const std::string &type_out, double *out );

		/**
		 * Convert an array of numbers from one unit to another.
		 * @param in 		A pointer to an array of double values
		 * @param nSamples 	The number of consecutive samples to convert
		 * @param type_in	String of the units to convert from
		 * @param type_out	The units to convert to
		 * @param out		A pointer to a preallocated array to hold the converted values
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 * @see units_t
		 */
		static void convert( const double *in, unsigned int nSamples,
			const std::string &type_in, units_t type_out, double *out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	String of the units to convert from
		 * @param type_out	The units to convert to
		 * @return 			The converted value
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 */
		static double convert( double in, const std::string &type_in, units_t type_out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	The units to convert from
		 * @param type_out	String of the units to convert to
		 * @return 			The converted value
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 */
		static double convert( double in, units_t type_in, const std::string &type_out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return 			The converted value
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 */
		static double convert( double in, units_t type_in, units_t type_out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	String of the units to convert from
		 * @param type_out	String of the units to convert to
		 * @return 			The converted value
		 * @throws invalid_conversion	if an undefined conversion is requested.
		 */
		static double convert( double in, const std::string &type_in, const std::string &type_out );

		/**
		 * Returns the string identification of the units type.
		 *
		 * @param type	The units constant to translate
		 * @return		The string identifying the constant
		 * @throws invalid_conversion if the constant is not recognized
		 */
		static std::string toString( units_t type );

		/**
		 * Returns the abbreviated string identification of the units type.
		 *
		 * @param type	The units constant to translate
		 * @return		The abbreviation identifying the constant
		 * @throws invalid_conversion if the constant is not recognized
		 */
		static std::string toStr( units_t type );

		/**
		 * Returns the units_t enum value associated with the supplied string.
		 *
		 * @param s 	The string to attempt to parse
		 * @return 		The enum value associated with the string
		 * @throws invalid_conversion if the string is not recognized
		 */
		static units_t fromString( std::string s = "" );

		/**
		 * Returns the units_t enum value associated with null units.
		 *
		 * @return 		The null enum value
		 */
		static units_t null();

		/**
		 * Prints a list of the recognized strings that can be translated to units_t
		 * values.
		 *
		 * @param o 	The ostream object to print the list to
		 */
		static void list_recognized_strings( std::ostream& o = std::cout );

		/**
		 * Retrieves a vector of recognized unit strings.
		 *
		 * @return 		A std::vector<std::string> list of recognized unit names
		 */
		static std::vector<std::string> get_recognized_strings();

	protected:

		static conversion_map_t map_;
		static string_to_units_map_t string_to_enum_map_;
		static units_to_string_map_t enum_to_string_map_, enum_to_abbr_map_;

		/**
		 * Generates a pair object associating the two units to be converted.
		 * @param t1		The unit to be converted from
		 * @param t2		The unit to be converted to
		 * @return 		A pair object that can be used as a map key
		 */
		static conversion_pair get_unit_pair_( units_t t1, units_t t2 );

		static void initialize_();
		static bool ready_();
	};
}



#endif
