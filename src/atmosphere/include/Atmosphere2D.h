#ifndef NCPAPROP_ATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE2D_H_INCLUDED

#include "AtmosphericModel.h"
#include "Atmosphere1D.h"
#include "NCPAUnits.h"
#include "NCPACommon.h"
#include <vector>
#include <climits>
#include <utility>

/**
 * Defines the interface for a 2D atmosphere (r,z).
 */

namespace NCPA { class Atmosphere2D; }
void swap( NCPA::Atmosphere2D &a, NCPA::Atmosphere2D &b );

namespace NCPA {

//	bool sort_profiles_by_range_( Atmosphere1D *p1, Atmosphere1D *p2 );

	class Atmosphere2D : public AtmosphericModel {

	public:
		Atmosphere2D();
		Atmosphere2D( Atmosphere2D &&in );
		virtual ~Atmosphere2D();
		friend void ::swap(Atmosphere2D &a, Atmosphere2D &b);

		// These will call the pure virtual method with type conversions
		// as appropriate
		virtual void set_range_units( const std::string &new_units );
		virtual void convert_range_units( const std::string &new_units );
		virtual void set_altitude_units( const std::string &s );
		virtual void convert_altitude_units( const std::string &new_units );
		virtual void set_property_units( const std::string &key,
				const std::string &new_units );
		virtual void convert_property_units( const std::string &key,
				const std::string &new_units );
		virtual void set_maximum_valid_range( double maxrange,
				const std::string &units = "" );
		virtual void set_maximum_valid_range( double maxrange,
				units_t units = NULL_UNIT );
		virtual void print_atmosphere(
				const std::vector< std::string >& columnorder,
				double range = 0.0, const std::string &altitude_key = "Z",
				std::ostream& os = std::cout );

		// backward compatibility methods
		virtual double get_first_derivative( double range, const std::string &key,
				double altitude );
		virtual double get_second_derivative( double range, const std::string &key,
				double altitude );

		// interface: pure virtual methods
		virtual double get( const std::string &key, double r, double z ) 		= 0;
		virtual double get_as( const std::string key, double r, double z,
				units_t units) 													= 0;
		virtual double get_as( const std::string key, double r, double z,
				const std::string &units) 										= 0;
		virtual double get_ddr( const std::string key, double r, double z ) 	= 0;
		virtual double get_ddz( const std::string key, double r, double z ) 	= 0;
		virtual double get_d2dr2( const std::string key, double r, double z ) 	= 0;
		virtual double get_d2dz2( const std::string key, double r, double z ) 	= 0;
		virtual double get_d2drdz( const std::string key, double r, double z ) 	= 0;

		virtual bool contains( const std::string &key ) const 					= 0;
		virtual bool contains_key( const std::string &key ) const 				= 0;
		virtual bool contains_scalar( const std::string &key ) const 			= 0;
		virtual bool contains_vector( const std::string &key ) const 			= 0;

		virtual std::vector<std::string> get_keys() const						= 0;
		virtual std::vector<std::string> get_vector_keys() const				= 0;
		virtual std::vector<std::string> get_scalar_keys() const				= 0;

		virtual void copy_property( const std::string &old_key,
				const std::string &new_key ) 									= 0;
		virtual void remove_property( const std::string &key ) 					= 0;

		virtual units_t get_range_units() const 								= 0;
		virtual void set_range_units( units_t new_units ) 						= 0;
		virtual void convert_range_units( units_t new_units ) 					= 0;

		virtual units_t get_altitude_units( double range ) const 				= 0;
		virtual void set_altitude_units( units_t new_units ) 					= 0;
		virtual void convert_altitude_units( units_t new_units ) 				= 0;

		virtual units_t get_property_units( const std::string &key ) const 		= 0;
		virtual void set_property_units( units_t new_units )					= 0;
		virtual void convert_property_units( const std::string &key,
				units_t new_units ) 											= 0;

		virtual ScalarWithUnits get_minimum_altitude( double range ) const		= 0;
		virtual ScalarWithUnits get_maximum_altitude( double range ) const 		= 0;

		virtual void get_minimum_altitude_limits( ScalarWithUnits &lowlimit,
				ScalarWithUnits &highlimit ) const 								= 0;
		virtual void get_maximum_altitude_limits( ScalarWithUnits &lowlimit,
				ScalarWithUnits &highlimit ) const 								= 0;
		virtual ScalarWithUnits get_maximum_valid_range() const 				= 0;
		virtual void set_maximum_valid_range( const ScalarWithUnits &s )		= 0;











//
//
//
//		virtual double get_interpolated_ground_elevation( double range );
//
//		virtual double get_interpolated_ground_elevation_first_derivative( double range );
//
//		virtual double get_interpolated_ground_elevation_second_derivative( double range );











//		size_t get_profile_index( double range, units_t units = NULL_UNIT );
//		size_t get_profile_index( double range, const std::string &units );

//		void get_property_vector( double range, const std::string &key, double *buffer,
//			units_t *buffer_units );
//		void get_property_vector( double range, const std::string &key, double *buffer );

//		units_t get_range_units() const;






//		void insert_profile( const Atmosphere1D *profile, double range,
//				const std::string &units = "" );
//		void insert_profile( const Atmosphere1D *profile, double range,
//				units_t units = NULL_UNIT );

//		std::vector< NCPA::Atmosphere1D * >::iterator last_profile();

//		size_t nz( double range, units_t units = NULL_UNIT );
//		size_t nz( double range, const std::string &units );



//		void read_elevation_from_file( const std::string &filename );






//		void set_insert_range_units( units_t u );
//		void set_insert_range_units( const std::string &u );



//		void sort_profiles();


	// bulk calculations
//		void calculate_attenuation( const std::string &new_key,
//			const std::string &temperature_key, const std::string &pressure_key,
//			const std::string &density_key, double freq, double tweak_factor = 1.0 );
//		void calculate_attenuation( const std::string &new_key,
//			const std::string &temperature_key, const std::string &pressure_key,
//			const std::string &density_key, const std::string &humidity_key,
//			double freq, double tweak_factor = 1.0 );
//		void read_attenuation_from_file( const std::string &new_key,
//			const std::string &filename );    // apply universally
		
//	protected:
//		void clear_last_index_();
//		void set_last_index_( size_t ind );
//		void calculate_midpoints_();
//		void generate_ground_elevation_spline_();
//		void free_ground_elevation_spline_();
//		void setup_ground_elevation_spline_from_profiles_();
//		void setup_ground_elevation_spline_from_file_();
//		units_t check_and_override_range_units_( units_t new_units );



		// for keeping track of the last index calculated, speeds things up a little
//		size_t last_index_ = 0;
//		units_t range_units_ = NULL_UNIT;
//		bool sorted_ = false;

//		ScalarWithUnits last_index_min_range_, last_index_max_range_;

		// data storage
//		std::vector< Atmosphere1D * > profiles_;
//		std::vector< ScalarWithUnits > midpoints_;
//
//		// ground elevation splines
//		double *topo_ground_heights_ = nullptr;
//		double *topo_ranges_ = nullptr;
//		gsl_interp_accel *topo_accel_ = nullptr;
//		gsl_spline *topo_spline_ = nullptr;
//		bool override_profile_z0_ = false;
//		std::string ground_elevation_file_ = "";
//
//		ScalarWithUnits max_valid_range_;

	};


}







#endif
