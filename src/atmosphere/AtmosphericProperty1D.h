#ifndef NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED
#define NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED

//#include <map>
//#include <stack>
#include <utility>
#include "NCPAUnits.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#if GSL_MAJOR_VERSION > 1
#define ATMOSPHERIC_INTERPOLATION_TYPE gsl_interp_steffen
#else
#define ATMOSPHERIC_INTERPOLATION_TYPE gsl_interp_cspline
#endif

namespace NCPA {

	class AtmosphericProperty1D : public std::pair< NCPA::VectorWithUnits, NCPA::VectorWithUnits > {
	protected:
//		double *z_;
		gsl_interp_accel *accel_ = NULL;
		gsl_spline *spline_ = NULL;
//		units_t z_units_;

//		NCPA::VectorWithUnits z_;

		int check_altitude_( double z_req ) const;
		void build_splines_();
		void delete_splines_();

	public:
		AtmosphericProperty1D();
		AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units );
		AtmosphericProperty1D( const AtmosphericProperty1D &source );
		virtual ~AtmosphericProperty1D();

		void convert_altitude_units( units_t new_units );
		units_t get_altitude_units() const;
		units_t get_units() const;
		virtual void convert_units( units_t new_units );

		void resample( double new_dz );

		void get_altitude_vector( double *buffer, units_t *buffer_units ) const;
		void get_altitude_vector( double *buffer ) const;
		void get_vector( double *buffer, units_t *buffer_units ) const;
		void get_vector( double *buffer ) const;
		double get( double altitude ) const;
		double get_first_derivative( double altitude ) const;
		double get_second_derivative( double altitude ) const;
	};

}

#endif
