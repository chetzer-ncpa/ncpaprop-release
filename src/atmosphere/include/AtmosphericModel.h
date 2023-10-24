#ifndef NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED
#define NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED

#include <string>
#include "NCPACommon.h"
#include "NCPAUnits.h"

#ifndef NCPA__GAMMA_FOR_C
#define NCPA__GAMMA_FOR_C 1.4
#endif

#ifndef NCPA__R_FOR_C
#define NCPA__R_FOR_C 287.0
#endif

namespace NCPA { class AtmosphericModel; }
void swap( NCPA::AtmosphericModel &a, NCPA::AtmosphericModel &b ) noexcept;

namespace NCPA {

	typedef enum atmospheric_stability_t : unsigned int {
		ATMOSPHERE_STABILITY_INVALID = 0,
		ATMOSPHERE_STABILITY_UNSTABLE,
		ATMOSPHERE_STABILITY_NEUTRAL,
		ATMOSPHERE_STABILITY_STABLE
	} atmospheric_stability_t;

	class AtmosphericModel {
	public:

		// swap function, does nothing cause there's (currently) nothing to swap
		friend void ::swap( AtmosphericModel &a, AtmosphericModel &b ) noexcept;
		virtual ~AtmosphericModel() = default;


		// t in K, returns in m/s
		static double soundspeed_from_temperature( double t,
				units_t t_units = Units::fromString("K"),
				units_t c_units = Units::fromString("m/s") );
		static double soundspeed_from_temperature( double t,
				const std::string &t_units,	const std::string &c_units );

		// c in m/s, returns in K
		static double temperature_from_soundspeed( double c,
				units_t c_units = Units::fromString("m/s"),
				units_t t_units = Units::fromString("K") );
		static double temperature_from_soundspeed( double c,
				const std::string &t_units,
				const std::string &c_units );

		// p in Pa, d in kg/m3, returns in m/s
//		static double soundspeed_from_pressure_density( double p, double d );
		static double soundspeed_from_pressure_density( double p, units_t p_units,
				double d, units_t d_units, units_t c_units );
		static double soundspeed_from_pressure_density( double p, const std::string &p_units,
				double d, const std::string &d_units, const std::string &c_units );


		static double density_from_temperature_pressure(
				double t, NCPA::units_t T_units,
				double p, NCPA::units_t P_units,
				NCPA::units_t D_units );
		static double density_from_temperature_pressure(
				double t, const std::string &T_units,
				double p, const std::string &P_units,
				const std::string &D_units );


		// z in km, t in K, p in Pa, d in kg/m3, freq in Hz
		static double attenuation_sutherland_bass(
			double z_km, double t_K, double p_Pa, double d_kgpm3, double freq );
		static double attenuation_from_temperature_pressure_density(
			double z_km, double t_K, double p_Pa, double d_kgpm3, double freq );
		static double attenuation_from_temperature_pressure_humidity(
			double z_km, double t_K, double p_Pa, double humidity_dec, double freq );
//		static double density_from_temperature_pressure( double t, double p );

		static atmospheric_stability_t get_stability_type( std::string s );

		virtual void calculate_sound_speed_from_temperature( const std::string &new_key,
			const std::string &temperature_key, NCPA::units_t speed_units ) = 0;
		virtual void calculate_sound_speed_from_pressure_and_density(
			const std::string &new_key,
			const std::string &pressure_key, const std::string &density_key,
			NCPA::units_t speed_units ) = 0;
		virtual void calculate_density_from_temperature_and_pressure(
			const std::string &new_key, const std::string &temperature_key,
			const std::string &pressure_key, units_t density_units ) = 0;
		virtual void calculate_effective_sound_speed( const std::string &new_key,
			const std::string &sound_speed_key, const std::string &wind_component_key ) = 0;
		virtual void calculate_wind_component( const std::string &new_key,
			const std::string &wind_speed_key, const std::string &wind_direction_key,
			double azimuth ) = 0;
		virtual void calculate_wind_direction( const std::string &new_key,
			const std::string &we_wind_speed_key, const std::string &sn_wind_speed_key,
			units_t direction_units = Units::fromString("Degrees clockwise from north")
			) = 0;
		virtual void calculate_wind_speed( const std::string &new_key,
			const std::string &we_wind_speed_key,
			const std::string &sn_wind_speed_key ) = 0;


	protected:
		std::string remove_underscores( const std::string &input ) const;
	};

}










#endif
