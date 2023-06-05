#ifndef NCPA_ATTENUATION_H_INCLUDED
#define NCPA_ATTENUATION_H_INCLUDED

#include <cstring>
#include <map>
#include <vector>

// Other packages
#include "NCPAUnits.h"

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_Z
#define NCPA_ATTENUATION_DEFAULT_UNITS_Z NCPA::UNITS_DISTANCE_KILOMETERS
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_T
#define NCPA_ATTENUATION_DEFAULT_UNITS_T NCPA::UNITS_TEMPERATURE_KELVIN
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_D
#define NCPA_ATTENUATION_DEFAULT_UNITS_D NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_P
#define NCPA_ATTENUATION_DEFAULT_UNITS_P NCPA::UNITS_PRESSURE_PASCALS
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_H
#define NCPA_ATTENUATION_DEFAULT_UNITS_H NCPA::UNITS_NONE
#endif

namespace NCPA {

	enum attenuation_model_parameter_label_t : unsigned int {
		ATTENUATION_MODEL_PARAMETER_ALTITUDE,		// default km
		ATTENUATION_MODEL_PARAMETER_TEMPERATURE,	// default K
		ATTENUATION_MODEL_PARAMETER_DENSITY,		// default kg/m3
		ATTENUATION_MODEL_PARAMETER_HUMIDITY,		// default pct
		ATTENUATION_MODEL_PARAMETER_PRESSURE		// default Pa
	};

	typedef NCPA::

	typedef std::pair<attenuation_model_parameter_label_t,NCPA::VectorWithUnits *> attenuation_parameter_pair_t;
	typedef std::map<NCPA::attenuation_model_parameter_label_t, NCPA::VectorWithUnits *> attenuation_parameter_map_t;

	class AttenuationModel {
		public:

			virtual ~AttenuationModel();

			virtual void set_vertical_scale(size_t nz, const double *z, NCPA::units_t units = NCPA_ATTENUATION_DEFAULT_UNITS_Z);
			virtual void set_vertical_scale(const NCPA::VectorWithUnits &z);

			virtual void set_vertical_parameter(attenuation_model_parameter_label_t param, size_t nz, const double *p,
					NCPA::units_t units);
			virtual void set_vertical_parameter(attenuation_model_parameter_label_t param,
					const NCPA::VectorWithUnits &p);

			// pure virtual methods: must be overridden
			virtual double attenuation(double f, double z, NCPA::units_t z_units = UNITS_NONE) = 0;
			virtual void attenuation(size_t n, double f, double *z, double *alpha, NCPA::units_t z_units = UNITS_NONE) = 0;
			virtual units_t get_calculation_units( attenuation_model_parameter_label_t param ) const = 0;
			virtual const std::vector<attenuation_model_parameter_label_t> get_required_parameters() const = 0;

		protected:
			attenuation_parameter_map_t parameter_map;
			void _set_parameter(attenuation_model_parameter_label_t param, NCPA::VectorWithUnits *p);
	};
}

#endif
