#ifndef NCPAPROP_PROFILEGRIDATMOSPHERE3D_H_INCLUDED
#define NCPAPROP_PROFILEGRIDATMOSPHERE3D_H_INCLUDED

#include "NCPACommon.h"
#include "Atmosphere3D.h"
#include <climits>
#include <map>



namespace NCPA {


	class ProfileGridAtmosphere3D : public Atmosphere3D {
		public:
			ProfileGridAtmosphere3D( const std::string &summary_file, size_t nz, double *zvec,
				size_t precision, const std::string &header_file );
			ProfileGridAtmosphere3D( const std::string &summary_file, size_t nz, double *zvec,
				size_t precision );
			ProfileGridAtmosphere3D( const std::string &summary_file, size_t nz, double *zvec,
				const std::string &header_file );
			ProfileGridAtmosphere3D( const std::string &summary_file, size_t nz, double *zvec );
			virtual ~ProfileGridAtmosphere3D();

			// setup and building
			void read_atmosphere_from_file( const std::string &summary_file,
				size_t nz, double *zvec, size_t precision, const std::string &header_file );
			virtual void convert_range_units( units_t new_units );
			virtual void convert_altitude_units( units_t new_units );
			virtual void convert_property_units( const std::string &key, units_t new_units );


			// querying
			virtual double get_maximum_valid_range( double azimuth ) const;
			virtual bool contains_property( const std::string &key ) const;
			virtual bool contains_vector( const std::string &key ) const;
			virtual bool contains_scalar( const std::string &key ) const;
			virtual units_t get_range_units() const;
			virtual units_t get_altitude_units() const;
			virtual units_t get_property_units( const std::string &key ) const;
			virtual double get_minimum_altitude( double x, double y ) const;
			virtual void get_minimum_altitude_limits(
				double &minlimit, double &maxlimit ) const;
			virtual void get_maximum_altitude_limits(
				double &minlimit, double &maxlimit ) const;
			// virtual void get_dimensions( size_t &nx, size_t &ny, size_t &nz ) const;

			// data retrieval
			virtual double get( double x, double y, const std::string &key ) const;   // get scalar quantity
			virtual double get( double x, double y, double z, const std::string &key ) const;  // get vector quantity
			virtual double get_derivative( double x, double y, const std::string &key,
				size_t order, deriv_t *directions ) const;
			virtual double get_derivative( double x, double y, double z,
				const std::string &key, size_t order, deriv_t *directions ) const;

			// assume that we're adding properties to an existing atmosphere, and these
			// are at the same grid points as the existing ones
			// @todo creation from scratch
			// manipulation
			virtual void add_property( const std::string &key, double ***prop,
				size_t nx, size_t ny, size_t nz, NCPA::units_t units );
			virtual void add_property( const std::string &key, double **prop,
				size_t nx, size_t ny, NCPA::units_t units );
			virtual NCPA::AtmosphericProperty3D *get_property( const std::string &key ) const;
			virtual void remove_property( const std::string &key );
			virtual void remove_vector_property( const std::string &key );
			virtual void remove_scalar_property( const std::string &key );
			virtual void get_property_template( const std::string &basis,
				size_t &nx, double *x, NCPA::units_t &x_units,
				size_t &ny, double *y, NCPA::units_t &y_units,
				double **&prop ) const;
			virtual void get_property_template( const std::string &basis,
				size_t &nx, double *x, NCPA::units_t &x_units,
				size_t &ny, double *y, NCPA::units_t &y_units,
				size_t &nz, double *z, NCPA::units_t &z_units,
				double ***&prop ) const;
			virtual void free_property_template(
				size_t nx, double *x,
				size_t ny, double *y,
				double **prop ) const;
			virtual void free_property_template(
				size_t nx, double *x,
				size_t ny, double *y,
				size_t nz, double *z,
				double ***prop ) const;
			virtual void copy_property( const std::string &old_key,
				const std::string &new_key );

			// bulk calculations
			virtual void calculate_sound_speed_from_temperature( const std::string &new_key,
				const std::string &temperature_key,  NCPA::units_t wind_units );
			virtual void calculate_sound_speed_from_pressure_and_density(
				const std::string &new_key, const std::string &pressure_key,
				const std::string &density_key, units_t wind_units );
			virtual void calculate_wind_speed( const std::string &new_key,
				const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key );
			virtual void calculate_wind_direction( const std::string &new_key,
				const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
				units_t direction_units = NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
			virtual void calculate_attenuation( const std::string &new_key,
				const std::string &temperature_key, const std::string &pressure_key,
				const std::string &density_key, double freq, double tweak_factor = 1.0 );
			virtual void calculate_wind_component( const std::string &new_key,
				const std::string &wind_speed_key, const std::string &wind_direction_key,
				double azimuth );
			virtual void calculate_effective_sound_speed( const std::string &new_key,
				const std::string &sound_speed_key, const std::string &wind_component_key );
			virtual void read_attenuation_from_file( const std::string &key,
				const std::string &filename );


		protected:
			std::map< std::string, VectorAtmosphericProperty3D *> vector_contents_;
			std::map< std::string, ScalarAtmosphericProperty3D *> scalar_contents_;

			void add_property( const std::string &key,
				NCPA::VectorAtmosphericProperty3D *prop );
			void add_property( const std::string &key,
				NCPA::ScalarAtmosphericProperty3D *prop );

	};

}





#endif   // NCPAPROP_ATMOSPHERE3D_H_INCLUDED
