/*
Classes, functions and constants for unit conversion.

Classes:
NCPA::UnitConverter: 	Handle all defined unit conversions.  Can convert single values or vectors.
			Will throw an exception if an undefined conversion is requested.

Constants:
NCPA::units_t	NCPA::UNITS_NONE
NCPA::units_t 	NCPA::UNITS_TEMPERATURE_CELSIUS
NCPA::units_t 	NCPA::UNITS_TEMPERATURE_KELVIN
NCPA::units_t 	NCPA::UNITS_TEMPERATURE_FAHRENHEIT
NCPA::units_t 	NCPA::UNITS_DISTANCE_METERS
NCPA::units_t 	NCPA::UNITS_DISTANCE_KILOMETERS
NCPA::units_t 	NCPA::UNITS_SPEED_METERS_PER_SECOND
NCPA::units_t 	NCPA::UNITS_SPEED_KILOMETERS_PER_SECOND
NCPA::units_t 	NCPA::UNITS_PRESSURE_PASCALS
NCPA::units_t 	NCPA::UNITS_PRESSURE_MILLIBARS
NCPA::units_t   NCPA::UNITS_PRESSURE_HECTOPASCALS
NCPA::units_t 	NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER
NCPA::units_t 	NCPA::UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER
NCPA::units_t	NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH
NCPA::units_t	NCPA::UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST
NCPA::units_t 	NCPA::UNITS_ANGLE_DEGREES
NCPA::units_t 	NCPA::UNITS_ANGLE_RADIANS

Examples:
	using namespace NCPA;

	// Convert a single value
	double temp_c = 50.0, temp_k;
	temp_k = Units::convert( temp_c, UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_KELVIN );
	// or
	Units::convert( &temp_c, 1, UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_KELVIN, &temp_k );

	// Convert a vector of values
	double Pa_vec[ 20 ] = { ... };
	double mbar_vec[ 20 ];
	Units::convert( Pa_vec, 20, UNITS_PRESSURE_PASCALS, UNITS_PRESSURE_MILLIBARS, mbar_vec );

	// Can also convert in-place
	double distance[ 300 ] = { ... };
	Units::convert( distance, 300, UNITS_DISTANCE_METERS, UNITS_DISTANCE_KILOMETERS, distance );


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

namespace NCPA {

	/**
	 * An enum of unit constants.
	 * Constants that can be used to identify or specify units.
	 */
	typedef enum units_t : unsigned int {
		UNITS_NONE = 0,						/**< Indicates no units */

		UNITS_TEMPERATURE_KELVIN,				/**< Temperature in Kelvin */
		UNITS_TEMPERATURE_CELSIUS,				/**< Temperature in Celsius */
		UNITS_TEMPERATURE_FAHRENHEIT,				/**< Temperature in Fahrenheit */

		UNITS_DISTANCE_METERS,					/**< Distance in meters */
		UNITS_DISTANCE_KILOMETERS,				/**< Distance in kilometers */

		UNITS_SPEED_METERS_PER_SECOND,				/**< Speed in m/s */
		UNITS_SPEED_KILOMETERS_PER_SECOND,			/**< Speed in km/s */

		UNITS_PRESSURE_PASCALS,					/**< Pressure in Pa */
		UNITS_PRESSURE_MILLIBARS,				/**< Pressure in mbar */
		UNITS_PRESSURE_HECTOPASCALS,			/**< Pressure in hPa */
		UNITS_PRESSURE_ATMOSPHERES,				/**< Pressure in atm */

		UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER,		/**< Density in kg/m^3 */
		UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER,		/**< Density in g/cm^3 */

		UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH,		/**< Direction in geographic azimuth */
		UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST,	/**< Direction in "math" convention */

		UNITS_ANGLE_DEGREES,					/**< Angles in degrees */
		UNITS_ANGLE_RADIANS						/**< Angles in radians */
	} units_t;
}


typedef std::pair< NCPA::units_t, NCPA::units_t > conversion_pair;
typedef double (*conversion_function)(double);
typedef std::map< conversion_pair, conversion_function > conversion_map_t;
typedef std::map< std::string, NCPA::units_t > string_to_units_map_t;
typedef std::map< NCPA::units_t, std::string > units_to_string_map_t;

namespace NCPA {

	/**
	 * A class for converting units.
	 */
	class Units {

	public:
		/**
		 * Convert an array of numbers from one unit to another.
		 * @param in 		A pointer to an array of double values
		 * @param nSamples 	The number of consecutive samples to convert
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @param out		A pointer to a preallocated array to hold the converted values
		 * @throws out_of_range	if an undefined conversion is requested.
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
		 * @throws out_of_range	if an undefined conversion is requested.
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
		 * @throws out_of_range	if an undefined conversion is requested.
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
		 * @throws out_of_range	if an undefined conversion is requested.
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
		 * @throws out_of_range	if an undefined conversion is requested.
		 */
		static double convert( double in, const std::string &type_in, units_t type_out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	The units to convert from
		 * @param type_out	String of the units to convert to
		 * @return 			The converted value
		 * @throws out_of_range	if an undefined conversion is requested.
		 */
		static double convert( double in, units_t type_in, const std::string &type_out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return 			The converted value
		 * @throws out_of_range	if an undefined conversion is requested.
		 */
		static double convert( double in, units_t type_in, units_t type_out );

		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	String of the units to convert from
		 * @param type_out	String of the units to convert to
		 * @return 			The converted value
		 * @throws out_of_range	if an undefined conversion is requested.
		 */
		static double convert( double in, const std::string &type_in, const std::string &type_out );

		/**
		 * Returns the string identification of the units type.
		 *
		 * @param type	The units constant to translate
		 * @return		The string identifying the constant
		 * @throws out_of_range if the constant is not recognized
		 */
		static std::string toString( units_t type );

		/**
		 * Returns the abbreviated string identification of the units type.
		 *
		 * @param type	The units constant to translate
		 * @return		The abbreviation identifying the constant
		 * @throws out_of_range if the constant is not recognized
		 */
		static std::string toStr( units_t type );

		/**
		 * Returns the units_t enum value associated with the supplied string.
		 *
		 * @param s 	The string to attempt to parse
		 * @return 		The enum value associated with the string
		 * @throws out_of_range if the string is not recognized
		 */
		static units_t fromString( std::string s );

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
