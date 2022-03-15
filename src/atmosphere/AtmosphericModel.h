#ifndef NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED
#define NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED

#include <string>
#include "ncpaprop_common.h"

namespace NCPA {

	typedef enum atmospheric_stability_t : unsigned int {
		ATMOSPHERE_STABILITY_INVALID = 0,
		ATMOSPHERE_STABILITY_UNSTABLE,
		ATMOSPHERE_STABILITY_NEUTRAL,
		ATMOSPHERE_STABILITY_STABLE
	} atmospheric_stability_t;

	class AtmosphericModel {
	public:
		// t in K, returns in m/s
		static double soundspeed_from_temperature( double t );

		// p in Pa, d in kg/m3, returns in m/s
		static double soundspeed_from_pressure_density( double p, double d );

		// z in km, t in K, p in Pa, d in kg/m3, freq in Hz
		static double attenuation_sutherland_bass(
			double z_km, double t_K, double p_Pa, double d_kgpm3, double freq );
		static double attenuation_from_temperature_pressure_density(
			double z_km, double t_K, double p_Pa, double d_kgpm3, double freq );
		static double attenuation_from_temperature_pressure_humidity(
			double z_km, double t_K, double p_Pa, double humidity_dec, double freq );

		static atmospheric_stability_t get_stability_type( std::string s );

		virtual ~AtmosphericModel() = default;
		virtual void calculate_sound_speed_from_temperature( const std::string &new_key,
			const std::string &temperature_key, NCPA::units_t speed_units ) = 0;
		virtual void calculate_sound_speed_from_pressure_and_density(
			const std::string &new_key,
			const std::string &pressure_key, const std::string &density_key,
			NCPA::units_t speed_units ) = 0;
	};

}










#endif
