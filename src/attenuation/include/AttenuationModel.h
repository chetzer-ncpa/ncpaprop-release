#ifndef NCPA_ATTENUATION_H_INCLUDED
#define NCPA_ATTENUATION_H_INCLUDED

#include <cstring>
#include <map>
#include <vector>

// Other packages
#include "NCPAUnits.h"
#include "NCPAInterpolation.h"

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_Z
#define NCPA_ATTENUATION_DEFAULT_UNITS_Z NCPA::units_t::DISTANCE_KILOMETERS
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_T
#define NCPA_ATTENUATION_DEFAULT_UNITS_T NCPA::units_t::TEMPERATURE_KELVIN
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_D
#define NCPA_ATTENUATION_DEFAULT_UNITS_D NCPA::units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_P
#define NCPA_ATTENUATION_DEFAULT_UNITS_P NCPA::units_t::PRESSURE_PASCALS
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_H
#define NCPA_ATTENUATION_DEFAULT_UNITS_H NCPA::units_t::RATIO_PERCENT
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_UNITS_A
#define NCPA_ATTENUATION_DEFAULT_UNITS_A NCPA::units_t::ATTENUATION_NEPERS_PER_METER
#endif

#ifndef NCPA_ATTENUATION_DEFAULT_INTERPOLATOR_TYPE
#define NCPA_ATTENUATION_DEFAULT_INTERPOLATOR_TYPE NCPA::interpolator1d_t::NCPA_1D_LINEAR
#endif

namespace NCPA { class AttenuationModel; }
void swap( NCPA::AttenuationModel &a, NCPA::AttenuationModel &b );

namespace NCPA {

	// typedefs and enums
	enum class attenuation_model_t : unsigned int {
		INVALID = 0,
		LOSSLESS,
		ISO9613_1,
		SUTHERLAND_BASS,
		TABULAR
	};

	enum class attenuation_parameter_t : unsigned int {
		ALTITUDE,		// default km
		TEMPERATURE,	// default K
		DENSITY,		// default kg/m3
		HUMIDITY,		// default pct
		PRESSURE,		// default Pa
		ATTENUATION		// default np/m
	};

	typedef std::map<attenuation_parameter_t,NCPA::VectorWithUnits>
		vector_attenuation_parameter_map_t;
	typedef std::map<attenuation_parameter_t,NCPA::ScalarWithUnits>
		scalar_attenuation_parameter_map_t;
	typedef std::map<attenuation_parameter_t,NCPA::Interpolator1D *>
		interpolator_map_t;

	class AttenuationModel {
		public:

			static AttenuationModel * build( attenuation_model_t t );
			static bool can_build( attenuation_model_t t );
			static std::string as_string( attenuation_model_t t );
			static std::string as_string( attenuation_parameter_t t );

			AttenuationModel() = default;
			AttenuationModel( const AttenuationModel &other );
			AttenuationModel( AttenuationModel &&other );
			friend void ::swap( AttenuationModel &a, AttenuationModel &b );
			virtual AttenuationModel* clone() const = 0;
			virtual ~AttenuationModel();

			virtual AttenuationModel* set_altitude_parameter(size_t nz, const double *z,
					NCPA::units_t units = NCPA_ATTENUATION_DEFAULT_UNITS_Z);
			virtual AttenuationModel* set_altitude_parameter(const std::vector<double> &z,
					NCPA::units_t units = NCPA_ATTENUATION_DEFAULT_UNITS_Z);
			virtual AttenuationModel* set_altitude_parameter(const NCPA::VectorWithUnits &z);
			virtual AttenuationModel* set_altitude_parameter(size_t nz, const double *z,
					const std::string &units);
			virtual AttenuationModel* set_altitude_parameter(const std::vector<double> &z,
					const std::string &units);


			virtual AttenuationModel* set_parameter(attenuation_parameter_t param, size_t nz,
					const double *p, NCPA::units_t units);
			virtual AttenuationModel* set_parameter(attenuation_parameter_t param, size_t nz,
					const double *p, const std::string &units);
			virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
					const std::vector<double> &p, NCPA::units_t units);
			virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
					const std::vector<double> &p, const std::string &units);
			virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
					const NCPA::VectorWithUnits &p);

			virtual AttenuationModel* set_parameter(attenuation_parameter_t param, double p,
					NCPA::units_t units);
			virtual AttenuationModel* set_parameter(attenuation_parameter_t param, double p,
					const std::string &units);
			virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
					const NCPA::ScalarWithUnits &p);

			virtual AttenuationModel* output_units( NCPA::units_t u );
			virtual NCPA::units_t output_units() const;

			virtual AttenuationModel* finalize();
			virtual AttenuationModel* interpolator( NCPA::interpolator1d_t inttype );
			virtual NCPA::interpolator1d_t interpolator() const;

			// pure virtual methods: must be overridden
			virtual double attenuation(double f, double z, NCPA::units_t z_units = units_t::NONE) = 0;
			virtual void attenuation(size_t n, double f, double *z, double *alpha, NCPA::units_t z_units = units_t::NONE) = 0;
			virtual units_t get_calculation_units( attenuation_parameter_t param ) const = 0;
			virtual const std::vector<attenuation_parameter_t> get_required_parameters() const = 0;
			virtual attenuation_model_t type() const = 0;

		protected:
			vector_attenuation_parameter_map_t vector_parameters_;
			scalar_attenuation_parameter_map_t scalar_parameters_;
			interpolator_map_t interpolators_;
			VectorWithUnits z_vector_;
			NCPA::interpolator1d_t interpolator_type_ = NCPA::interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR;
			NCPA::units_t output_units_ = NCPA::units_t::NONE,
						  native_output_units_ = NCPA::units_t::NONE;
			bool finalized_ = false;

			void set_parameter_(attenuation_parameter_t param,
					NCPA::VectorWithUnits p);
			void set_parameter_(attenuation_parameter_t param,
					NCPA::ScalarWithUnits p);
			void remove_parameter_( attenuation_parameter_t param );
			void clear_interpolators_();
	};
}

#endif
