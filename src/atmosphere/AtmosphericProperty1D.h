#pragma once
#include <map>
#include <stack>
#include "units.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#define NCPAPROP_ATMOSPHERICPROPERTY1D_DEFAULT_SAME_Z_TOLERANCE 1.0e-6
#define NCPAPROP_ATMOSPHERICPROPERTY1D_DEFAULT_SAME_Z_REPLACE false

namespace NCPA {

	class AtmosphericProperty1D : public VectorWithUnits {
		protected:
			double *z_ = nullptr;
			gsl_interp_accel *accel_ = nullptr;
			gsl_spline *spline_ = nullptr;
			units_t z_units_;

			int check_altitude_( double z_req ) const;
			void build_splines_();
			void delete_splines_();

		public:
			AtmosphericProperty1D();
			AtmosphericProperty1D( size_t n_points,
					double *altitude_points,
					units_t altitude_units,
					double *property_values,
					units_t property_units );
			AtmosphericProperty1D( const AtmosphericProperty1D &source );
			~AtmosphericProperty1D();

			double get( double altitude ) const;
			double get_first_derivative( double altitude ) const;
			double get_second_derivative( double altitude ) const;

			void convert_altitude_units( units_t new_units );
			units_t get_altitude_units() const;
			virtual void convert_units( units_t new_units );

			void resample( double new_dz );
			void add_point( double z,
					double f,
					bool replace =
						NCPAPROP_ATMOSPHERICPROPERTY1D_DEFAULT_SAME_Z_REPLACE,
					double tolerance =
						NCPAPROP_ATMOSPHERICPROPERTY1D_DEFAULT_SAME_Z_TOLERANCE );
			void add_points( size_t n_new_points,
					const double *z_points,
					const double *f_points,
					bool replace =
						NCPAPROP_ATMOSPHERICPROPERTY1D_DEFAULT_SAME_Z_REPLACE,
					double tolerance =
						NCPAPROP_ATMOSPHERICPROPERTY1D_DEFAULT_SAME_Z_TOLERANCE );

			void get_altitude_vector( double *buffer,
					units_t *buffer_units ) const;
	};

}
