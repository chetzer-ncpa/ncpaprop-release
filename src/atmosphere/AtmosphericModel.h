#ifndef NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED
#define NCPAPROP_ATMOSPHERICMODEL_H_INCLUDED

#include <string>
#include "ncpaprop_common.h"

namespace NCPA {

	class AtmosphericModel {
	public:
		// t in K, returns in m/s
		static double soundspeed_from_temperature( double t );

		// p in Pa, d in kg/m3, returns in m/s
		static double soundspeed_from_pressure_density( double p, double d );

		// z in km, t in K, p in Pa, d in kg/m3, freq in Hz
		static double attenuation_sutherland_bass(
			double z, double t, double p, double d, double freq );

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
