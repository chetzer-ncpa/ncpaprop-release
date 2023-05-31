#ifndef NCPA_ATTENUATION_H_INCLUDED
#define NCPA_ATTENUATION_H_INCLUDED

#include <cstring>
#include <map>
#include <vector>


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

	enum attenuation_model_parameter_t : unsigned int {
		ATTENUATION_MODEL_PARAMETER_ALTITUDE,		// default km
		ATTENUATION_MODEL_PARAMETER_TEMPERATURE,	// default K
		ATTENUATION_MODEL_PARAMETER_DENSITY,		// default kg/m3
		ATTENUATION_MODEL_PARAMETER_HUMIDITY,		// default pct
		ATTENUATION_MODEL_PARAMETER_PRESSURE		// default Pa
	};

	typedef std::pair<attenuation_model_parameter_t,NCPA::VectorWithUnits *> attenuation_parameter_pair_t;
	typedef std::map<NCPA::attenuation_model_parameter_t, NCPA::VectorWithUnits *> attenuation_parameter_map_t;

	class AttenuationModel {
		public:
			virtual ~AttenuationModel();

			void set_vertical_scale(size_t nz, const double *z, NCPA::units_t units = NCPA_ATTENUATION_DEFAULT_UNITS_Z);
			void set_vertical_scale(size_t nz, const NCPA::ScalarWithUnits *z);
			void set_vertical_scale(const NCPA::VectorWithUnits *z);

			void set_vertical_parameter(attenuation_model_parameter_t param, size_t nz, const double *p);
			void set_vertical_parameter(attenuation_model_parameter_t param, size_t nz, const double *p,
					NCPA::units_t units);
			void set_vertical_parameter(attenuation_model_parameter_t param, size_t nz,
					const NCPA::ScalarWithUnits *p);
			void set_vertical_parameter(attenuation_model_parameter_t param, const NCPA::VectorWithUnits *p);


			virtual double attenuation(double z, double f) = 0;
			virtual double attenuation(double z, double f,
					std::map<attenuation_model_parameter_t,double> params) = 0;
			virtual double attenuation(double z, NCPA::units_t z_units, double f) = 0;
			virtual double attenuation(NCPA::ScalarWithUnits *z, double f) = 0;
			virtual double attenuation(double z, NCPA::units_t z_units, double f,
					std::map<attenuation_model_parameter_t,double> params) = 0;
			virtual double attenuation(NCPA::ScalarWithUnits *z, double f,
					std::map<attenuation_model_parameter_t,double> params) = 0;

		protected:
			attenuation_parameter_map_t parameter_map;
			void _set_parameter(attenuation_model_parameter_t param, const NCPA::VectorWithUnits *p);







	};



}

#endif
