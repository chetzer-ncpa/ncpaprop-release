#ifndef NCPAPROP_ATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE2D_H_INCLUDED

#include "AtmosphericModel.h"
#include "Atmosphere1D.h"
#include <vector>
#include <climits>

namespace NCPA {

	bool sort_profiles_by_range_( Atmosphere1D *p1, Atmosphere1D *p2 );

	class Atmosphere2D : public AtmosphericModel {

	public:
		Atmosphere2D();
		virtual ~Atmosphere2D();

		// setup of profiles
		void insert_profile( const Atmosphere1D *profile, double range );
		void set_insert_range_units( units_t u );
		void sort_profiles();
		void convert_range_units( NCPA::units_t new_units );
		void set_maximum_valid_range( double maxrange );
		double get_maximum_valid_range() const;

		// profile manipulation
		void add_property( const std::string &key, size_t n_points, double *quantity_points,
			units_t quantity_units = NCPA::UNITS_NONE );
		void add_property( const std::string &key, double value,
			units_t units = NCPA::UNITS_NONE );
		void copy_vector_property( const std::string &old_key, const std::string &new_key );
		void copy_scalar_property( const std::string &old_key, const std::string &new_key );
		void remove_property( const std::string &key );


		// data retrieval, single values
		double get( double range, const std::string &key );    // retrieve scalar quantity from profile
		double get( double range, const std::string &key, double altitude );
		double get_first_derivative( double range, const std::string &key, double altitude );
		double get_second_derivative( double range, const std::string &key, double altitude );
		// @todo add radial derivatives


		// data retrieval, arrays
		size_t get_profile_index( double range );
		size_t nz( double range );
		void get_altitude_vector( double range, double *buffer, units_t *buffer_units );
		void get_property_vector( double range, const std::string &key, double *buffer,
			units_t *buffer_units );
		void get_altitude_vector( double range, double *buffer );
		void get_property_vector( double range, const std::string &key, double *buffer );
		units_t get_range_units() const;
		units_t get_altitude_units( double range );
		units_t get_property_units( double range, const std::string &key );
		bool contains_scalar( double range, const std::string &key );
		bool contains_vector( double range, const std::string &key );
		bool contains_key( double range, const std::string &key );
		virtual double get_interpolated_ground_elevation( double range );
		virtual double get_interpolated_ground_elevation_first_derivative( double range );
		virtual double get_interpolated_ground_elevation_second_derivative( double range );


		// metadata
		double get_minimum_altitude( double range );
		void get_minimum_altitude_limits( double &lowlimit, double &highlimit );
		double get_maximum_altitude( double range );
		void get_maximum_altitude_limits( double &lowlimit, double &highlimit );
		//double get_overall_maximum_altitude() const;

		// bulk calculations
		void calculate_density_from_temperature_and_pressure(
			const std::string &new_key, const std::string &temperature_key,
			const std::string &pressure_key, units_t density_units );
		void calculate_sound_speed_from_temperature( const std::string &new_key,
			const std::string &temperature_key, units_t wind_units );
		void calculate_sound_speed_from_pressure_and_density( const std::string &new_key,
			const std::string &pressure_key, const std::string &density_key,
			units_t wind_units );
		void calculate_wind_speed( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key );
		void calculate_wind_direction( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
			units_t direction_units = NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
		void calculate_attenuation( const std::string &new_key,
			const std::string &temperature_key, const std::string &pressure_key,
			const std::string &density_key, double freq, double tweak_factor = 1.0 );
		void calculate_attenuation( const std::string &new_key,
			const std::string &temperature_key, const std::string &pressure_key,
			const std::string &density_key, const std::string &humidity_key,
			double freq, double tweak_factor = 1.0 );
		void calculate_wind_component( const std::string &new_key,
			const std::string &wind_speed_key, const std::string &wind_direction_key,
			double azimuth );
		void calculate_effective_sound_speed( const std::string &new_key,
			const std::string &sound_speed_key, const std::string &wind_component_key );
		void convert_altitude_units( units_t new_units );
		void convert_property_units( const std::string &key, units_t new_units );
		void read_attenuation_from_file( const std::string &new_key,
			const std::string &filename );    // apply universally
		void read_elevation_from_file( const std::string &filename );
		void finalize_elevation_from_profiles();
		
		std::vector< NCPA::Atmosphere1D * >::iterator first_profile();
		std::vector< NCPA::Atmosphere1D * >::iterator last_profile();
	protected:
		void clear_last_index_();
		void set_last_index_( size_t ind );
		void calculate_midpoints_();
		void generate_ground_elevation_spline_();
		void free_ground_elevation_spline_();
		void setup_ground_elevation_spline_from_profiles_();
		void setup_ground_elevation_spline_from_file_();

		// data storage
		std::vector< Atmosphere1D * > profiles_;
		std::vector< double > midpoints_;

		// for keeping track of the last index calculated, speeds things up a little
		size_t last_index_;
		double last_index_min_range_, last_index_max_range_;
		units_t range_units_;
		bool sorted_;

		// ground elevation splines
		double *topo_ground_heights_;
		double *topo_ranges_;
		gsl_interp_accel *topo_accel_;
		gsl_spline *topo_spline_;
		bool override_profile_z0_ = false;
		std::string ground_elevation_file_;

		double max_valid_range_;

	};


}







#endif