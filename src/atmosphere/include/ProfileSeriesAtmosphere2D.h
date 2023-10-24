#ifndef NCPAPROP_PROFILESERIESATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_PROFILESERIESATMOSPHERE2D_H_INCLUDED

#include "Atmosphere2D.h"

namespace NCPA { class ProfileSeriesAtmosphere2D; }
void swap( NCPA::ProfileSeriesAtmosphere2D &a, NCPA::ProfileSeriesAtmosphere2D &b );

namespace NCPA {
	class ProfileSeriesAtmosphere2D : public Atmosphere2D {

	public:
		ProfileSeriesAtmosphere2D();
		ProfileSeriesAtmosphere2D( const std::string &filename, const std::string &headerfilename = "" );
		~ProfileSeriesAtmosphere2D();
		ProfileSeriesAtmosphere2D( const ProfileSeriesAtmosphere2D &other );
		ProfileSeriesAtmosphere2D( ProfileSeriesAtmosphere2D &&other );
		ProfileSeriesAtmosphere2D& operator=( ProfileSeriesAtmosphere other );

		void insert_profile( const Atmosphere1D *profile, double range,
				const std::string &units = "" );
		void insert_profile( const Atmosphere1D *profile, double range,
				units_t units = NULL_UNIT );

		void get_altitude_vector( double range, double *buffer, units_t *buffer_units );
		void get_altitude_vector( double range, double *buffer );


		// Atmosphere2D interface methods
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


	};
}




#endif
