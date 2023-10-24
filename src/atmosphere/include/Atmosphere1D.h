#ifndef NCPAPROP_ATMOSPHERE1D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE1D_H_INCLUDED

#include "NCPAInterpolation.h"
#include "AtmosphericModel.h"
#include "AtmosphericProperty1D.h"
#include "geographic.h"
#include "units.h"
#include <string>
#include <map>
#include <stack>
#include <fstream>
#include <algorithm>

#ifndef NCPA_ATMOSPHERE_SCIENTIFIC_PRECISION
#define NCPA_ATMOSPHERE_SCIENTIFIC_PRECISION 3
#endif

#ifndef NCPA_ATMOSPHERE_FIXED_PRECISION
#define NCPA_ATMOSPHERE_FIXED_PRECISION 3
#endif

#ifndef NCPA_ATMOSPHERE_FIXED_MAX_LOG
#define NCPA_ATMOSPHERE_FIXED_MAX_LOG 4
#endif

#ifndef NCPA_ATMOSPHERE_FIELD_WIDTH
#define NCPA_ATMOSPHERE_FIELD_WIDTH std::max((NCPA_ATMOSPHERE_SCIENTIFIC_PRECISION+8),(NCPA_ATMOSPHERE_FIXED_MAX_LOG+NCPA_ATMOSPHERE_FIXED_PRECISION+2))
#endif

namespace NCPA { class Atmosphere1D; }

void swap( NCPA::Atmosphere1D&, NCPA::Atmosphere1D& ) noexcept;

typedef std::map< std::string, NCPA::AtmosphericProperty1D * > vector_atmospheric_property_map_t;
typedef std::map< std::string, NCPA::ScalarWithUnits * > scalar_atmospheric_property_map_t;

namespace NCPA {

	/**
	A 1-dimensional atmospheric profile.
	A class that represents a 1-dimensional atmospheric model with properties that are either scalars
	or vectors that depend on altitude.  This allows for storage of, for example, a reference 
	ground altitude that allows zero to be used for sea level and, therefore, multiple profiles
	to be compared to one another using a common datum.  All properties have associated units
	that can be tracked and converted internally.
	*/
	class Atmosphere1D : public AtmosphericModel {

	public:

		Atmosphere1D();

		/**
		Constructs an Atmosphere1D object using a vector of altitudes, to which corresponding
		properties can then be added.
		@brief Constructor from altitude vector.
		@param n_altitude_points The number of points in the altitude vector
		@param altitude_points A pointer to the first altitude point
		@param altitude_units A Units object indicating the distance units used.
		*/
		Atmosphere1D( size_t n_altitude_points, double *altitude_points,
			units_t altitude_units );

		/**
		Constructs an Atmosphere1D object using a vector of altitudes, to which corresponding
		properties can then be added.
		@brief Constructor from altitude vector.
		@param n_altitude_points The number of points in the altitude vector
		@param altitude_points A pointer to the first altitude point
		@param altitude_units A Units object indicating the distance units used.
		*/
		Atmosphere1D( size_t n_altitude_points, double *altitude_points,
			const std::string &altitude_units );

		/**
		Constructs an Atmosphere1D object using formatted input from a file.
		@brief Constructor from formatted file.
		@param filename The name of the file to be read from.
		*/
		Atmosphere1D( const std::string &filename, const std::string &headerfilename = "" );

		/**
		Copy constructor.
		@brief Copy constructor.
		@param source The Atmosphere1D object to be copied.
		*/
		Atmosphere1D( const Atmosphere1D &source );

		/**
		Move constructor.
		@brief Move constructor.
		@param source The object to be moved.
		 */
		Atmosphere1D( Atmosphere1D &&source ) noexcept;

		/**
		Cleans up all dynamically-allocated memory and calls its superclass destructor.
		@brief Destructor
		*/
		~Atmosphere1D();

		/**
		 * Swap function.
		 * @brief Swaps contents of two Atmosphere1D objects A and B.
		 * @param a Object A
		 * @param b Object B
		 */
		friend void ::swap(Atmosphere1D &a, Atmosphere1D &b) noexcept;

		/**
		 * Assignment operator.
		 * @param a The object to be assigned to this.
		 */
		Atmosphere1D &operator=(Atmosphere1D a);

		/**
		Initializes the altitude (independent) vector.  Calls cleanup()
		first.
		*/
		void reset_altitude_vector( size_t n_altitude_points,
			double *altitude_points, units_t altitude_units );

		/**
		Initializes the altitude (independent) vector.  Calls cleanup()
		first.
		*/
		void reset_altitude_vector( size_t n_altitude_points,
			double *altitude_points, const std::string &altitude_units );

		/**
		Initializes the altitude (independent) vector.  Calls cleanup()
		first.
		*/
		void reset_altitude_vector( const NCPA::VectorWithUnits &v );



		/**
		Reads a formatted profile from an input stream, normally a std::ifstream representing a file.
		@brief Read a formatted atmospheric profile from a stream.
		@param in The input stream to read the profile from.  Should already be open.
		*/
		void read_values_from_stream( std::istream& in );

		/**
		Reads a formatted header from an input stream
		*/
		void read_header_from_stream( std::istream& in );

		/**
		Adds a property, including units, that depends on the altitude to the atmosphere.
		@brief Adds a vector property to the atmosphere.
		@param key The name by which the property will be referred.  Needs to be unique.
		@param n_points The number of points in the property vector.
		@param quantity_points A pointer to the first value in the property vector.
		@param quantity_units The units of the quantity.
		@throws std::runtime_error if the key name is already in use
		@throws std::runtime_error if the number of points does not match the altitude vector
		*/
		void add_property( const std::string &key, size_t n_points, double *quantity_points,
			units_t quantity_units = NCPA::units_t::NONE );

		/**
		Adds a property, including units, that depends on the altitude to the atmosphere.
		@brief Adds a vector property to the atmosphere.
		@param key The name by which the property will be referred.  Needs to be unique.
		@param n_points The number of points in the property vector.
		@param quantity_points A pointer to the first value in the property vector.
		@param quantity_units The units of the quantity.
		@throws std::runtime_error if the key name is already in use
		@throws std::runtime_error if the number of points does not match the altitude vector
		*/
		void add_property( const std::string &key, size_t n_points, double *quantity_points,
			const std::string &quantity_units = "" );

		/**
		Adds a property, including units, that depends on the altitude to the atmosphere.
		@brief Adds a vector property to the atmosphere.
		@param key The name by which the property will be referred.  Needs to be unique.
		@param p The Vector to insert
		@throws std::runtime_error if the key name is already in use
		@throws std::runtime_error if the number of points does not match the altitude vector
		*/
		void add_property( const std::string &key, const NCPA::VectorWithUnits &v );

		/**
		Adds a scalar property, including units, to the atmosphere.
		@brief Adds a scalar property to the atmosphere.
		@param key The name by which the property will be referred.  Needs to be unique.
		@param value The value of the property.
		@param units The units of the property.
		@throws std::runtime_error if the key name is already in use
		*/
		void add_property( const std::string &key, double value,
			units_t units = NCPA::units_t::NONE );

		/**
		Adds a scalar property, including units, to the atmosphere.
		@brief Adds a scalar property to the atmosphere.
		@param key The name by which the property will be referred.  Needs to be unique.
		@param s The Scalar object to add (copied)
		@throws std::runtime_error if the key name is already in use
		*/
		void add_property( const std::string &key, const NCPA::ScalarWithUnits &s );

		/**
		Adds a scalar property, including units, to the atmosphere.
		@brief Adds a scalar property to the atmosphere.
		@param key The name by which the property will be referred.  Needs to be unique.
		@param value The value of the property.
		@param units The units of the property.
		@throws std::runtime_error if the key name is already in use
		*/
		void add_property( const std::string &key, double value,
			const std::string &units = "" );


		void copy_vector_property( const std::string &old_key, const std::string &new_key );
		void copy_scalar_property( const std::string &old_key, const std::string &new_key );

		NCPA::AtmosphericProperty1D *get_vector_property_object( const std::string &key );
		const NCPA::AtmosphericProperty1D * const get_const_vector_property_object( const std::string &key ) const;
		NCPA::ScalarWithUnits *get_scalar_property_object( const std::string &key );
		const NCPA::ScalarWithUnits * const get_const_scalar_property_object( const std::string &key ) const;

		/**
		Removes a scalar or vector property from the atmosphere.  Has no effect if the property does not exist.
		@brief Removes a property.
		@param key The key for the property to be removed.
		*/
		void remove_property( const std::string &key );
		void reset_splines();


		double get( const std::string &key ) const;    // scalars
		double get( const std::string &key, double altitude ) const;
		double get_first_derivative( const std::string &key, double altitude ) const;
		double get_second_derivative( const std::string &key, double altitude ) const;

		double get_as( const std::string &key, units_t u ) const;    // scalars
		double get_as( const std::string &key, double altitude, units_t u ) const;
		double get_first_derivative_as( const std::string &key, double altitude,
				units_t u ) const;
		double get_second_derivative_as( const std::string &key, double altitude,
				units_t u ) const;

		double get_as( const std::string &key, const std::string &u ) const;    // scalars
		double get_as( const std::string &key, double altitude, const std::string &u ) const;
		double get_first_derivative_as( const std::string &key, double altitude,
				const std::string &u ) const;
		double get_second_derivative_as( const std::string &key, double altitude,
				const std::string &u ) const;


		size_t get_basis_length() const;
		size_t nz() const;
		void get_altitude_vector( double *buffer, units_t &buffer_units ) const;
		void get_property_vector( const std::string &key, double *buffer,
			units_t &buffer_units ) const;
		void get_altitude_vector( double *buffer ) const;
		void get_property_vector( const std::string &key, double *buffer ) const;
		units_t get_altitude_units() const;
		units_t get_property_units( const std::string &key ) const;


		/**
		Returns the minimum valid altitude in the atmospheric model, i.e. the first value in the altitude
		vector, in the current units.
		@brief Returns the minimum valid altitude.
		@returns The minimum valid altitude in the current units.
		*/
		double get_minimum_altitude() const;
		double get_maximum_altitude() const;

		// derived quantities
		void calculate_sound_speed_from_temperature( const std::string &new_key,
			const std::string &temperature_key, units_t speed_units );
		void calculate_sound_speed_from_temperature( const std::string &new_key,
			const std::string &temperature_key, const std::string &speed_units );

		void calculate_temperature_from_sound_speed( const std::string &new_key,
			const std::string &speed_key, units_t temp_units );
		void calculate_temperature_from_sound_speed( const std::string &new_key,
			const std::string &speed_key, const std::string &temp_units );

		void calculate_sound_speed_from_pressure_and_density( const std::string &new_key,
			const std::string &pressure_key, const std::string &density_key,
			units_t speed_units );
		void calculate_sound_speed_from_pressure_and_density( const std::string &new_key,
			const std::string &pressure_key, const std::string &density_key,
			const std::string &speed_units );

		void calculate_density_from_temperature_and_pressure(
			const std::string &new_key, const std::string &temperature_key,
			const std::string &pressure_key, units_t density_units );
		void calculate_density_from_temperature_and_pressure(
			const std::string &new_key, const std::string &temperature_key,
			const std::string &pressure_key, const std::string &density_units );

		void calculate_wind_speed( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key );
		void calculate_wind_direction( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
			units_t direction_units = NCPA::units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
		void calculate_wind_component( const std::string &new_key,
			const std::string &wind_speed_key, const std::string &wind_direction_key,
			double azimuth );
		void calculate_effective_sound_speed( const std::string &new_key,
			const std::string &sound_speed_key, const std::string &wind_component_key );

		// old version, no humidity option
		void calculate_attenuation( const std::string &new_key,
			const std::string &temperature_key, const std::string &pressure_key,
			const std::string &density_key, double freq, double tweak_factor = 1.0 );
		void calculate_attenuation( const std::string &new_key,
			const std::string &temperature_key, const std::string &pressure_key,
			const std::string &density_key, const std::string &humidity_key,
			double freq, double tweak_factor = 1.0 );
		void read_attenuation_from_file( const std::string &new_key,
			const std::string &filename );
		
		void convert_altitude_units( units_t new_units );
		void convert_property_units( const std::string &key, units_t new_units );
		//units_t get_property_units( std::string key );

		void resample( double new_dz, Interpolator1D *interp );
		void resample( VectorWithUnits &new_z, Interpolator1D *interp );

		std::vector< std::string > get_keys() const;
		std::vector< std::string > get_vector_keys() const;
		std::vector< std::string > get_scalar_keys() const;
		bool contains_scalar( const std::string &key ) const;
		bool contains_vector( const std::string &key ) const;
		bool contains_key( const std::string &key ) const;
		bool contains( const std::string &key ) const;

		void print_atmosphere(
				VectorWithUnits &z_values,
				const std::vector< std::string >& properties,
				std::ostream& os = std::cout,
				const std::string &z_label = "Z" ) const;
		void print_atmosphere(
				const std::vector< std::string >& properties,
				const std::string &altitude_key = "Z",
				std::ostream& os = std::cout ) const;
		std::string make_header_line( const std::string &key,
			size_t columnnumber ) const;
		std::string make_vector_header_line(
			const std::string &key, size_t columnnumber ) const;
		std::string make_scalar_header_line(
			const std::string &key ) const;

		static std::string format_header_line( size_t column,
				const std::string &key, const std::string &units,
				double val );
		static std::string format_header_line( size_t column,
				const std::string &key, const std::string &units,
				int val );
		static std::string format_header_line( size_t column,
				const std::string &key, const std::string &units );

		static std::string format_for_stream( double val );
		static std::string format_for_stream( int val );

	protected:
		// internal storage
		vector_atmospheric_property_map_t contents_;
		scalar_atmospheric_property_map_t scalar_contents_;
		
		NCPA::VectorWithUnits z_;
		std::vector< std::string > headerlines_;

		void cleanup_();
		void cleanup_vectors_();
		void cleanup_scalars_();
		void assert_key_does_not_exist_( const std::string &key ) const;
		void assert_key_exists_( const std::string &key ) const;
		void assert_vector_key_exists_( const std::string &key ) const;
		void assert_scalar_key_exists_( const std::string &key ) const;
		void do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits );

		// @todo: location
	};

}

#endif
