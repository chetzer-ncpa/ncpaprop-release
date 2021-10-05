#ifndef NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED
#define NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED

#include <map>
#include <stack>
#include "units.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"


namespace NCPA {

	class AtmosphericProperty1D : public VectorWithUnits {
	protected:
		double *z_;
		//NCPA::VectorWithUnits *z_vector_;
		gsl_interp_accel *accel_ = NULL;
		gsl_spline *spline_ = NULL;
		//std::stack< NCPA::units_t > z_units_;
		units_t z_units_;

		int check_altitude_( double z_req ) const;
		void build_splines_();
		void delete_splines_();

	public:
		AtmosphericProperty1D();
		AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units );
		AtmosphericProperty1D( const AtmosphericProperty1D &source );
		~AtmosphericProperty1D();

		double get( double altitude ) const;
		double get_first_derivative( double altitude ) const;
		double get_second_derivative( double altitude ) const;

		void convert_altitude_units( units_t new_units );
		units_t get_altitude_units() const;
		//void revert_altitude_units();
		virtual void convert_units( units_t new_units );
		//void revert_units();

		void resample( double new_dz );

		void get_altitude_vector( double *buffer, units_t *buffer_units ) const;
	};

}

#endif