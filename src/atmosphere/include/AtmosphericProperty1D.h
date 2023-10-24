#ifndef NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED
#define NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED

//#include <map>
//#include <stack>
#include <utility>
#include "NCPAUnits.h"
#include "NCPAInterpolation.h"

#ifndef NCPAPROP_DEFAULT_ATMOSPHERIC_PROPERTY_INTERPOLATION_TYPE
#define NCPAPROP_DEFAULT_ATMOSPHERIC_PROPERTY_INTERPOLATION_TYPE NCPA::interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR
#endif

namespace NCPA { class AtmosphericProperty1D; }
void swap( NCPA::AtmosphericProperty1D&, NCPA::AtmosphericProperty1D& ) noexcept;


namespace NCPA {
	typedef std::pair< NCPA::VectorWithUnits, NCPA::VectorWithUnits > vectorpair_t;

	class AtmosphericProperty1D : public vectorpair_t {
	public:
		AtmosphericProperty1D();
		AtmosphericProperty1D( size_t n_points,
				double *altitude_points, units_t altitude_units,
				double *property_values, units_t property_units,
				NCPA::interpolator1d_t interptype = NCPAPROP_DEFAULT_ATMOSPHERIC_PROPERTY_INTERPOLATION_TYPE);
		AtmosphericProperty1D( size_t n_points,
			double *altitude_points, const std::string &altitude_units,
			double *property_values, const std::string &property_units,
			NCPA::interpolator1d_t interptype = NCPAPROP_DEFAULT_ATMOSPHERIC_PROPERTY_INTERPOLATION_TYPE);
		AtmosphericProperty1D( NCPA::VectorWithUnits zvector, NCPA::VectorWithUnits propvector,
				NCPA::interpolator1d_t interptype = NCPAPROP_DEFAULT_ATMOSPHERIC_PROPERTY_INTERPOLATION_TYPE );
		AtmosphericProperty1D( const AtmosphericProperty1D &source );
		AtmosphericProperty1D( AtmosphericProperty1D &&source );
		virtual ~AtmosphericProperty1D();

		friend void ::swap( AtmosphericProperty1D &a, AtmosphericProperty1D &b ) noexcept;
		AtmosphericProperty1D& operator=( AtmosphericProperty1D source );

		virtual void convert_altitude_units( units_t new_units );
		virtual void convert_altitude_units( const std::string &new_units );
		virtual units_t get_altitude_units() const;
		virtual units_t get_units() const;
		virtual void convert_units( units_t new_units );
		virtual void convert_units( const std::string &new_units );
		virtual size_t size() const;

		virtual void reset_splines();

		virtual void get_altitude_vector( double *buffer, units_t &buffer_units );
		virtual void get_altitude_vector( double *buffer );
		virtual NCPA::VectorWithUnits get_altitude_vector();
		virtual void get_altitude_vector_as( double *buffer, units_t new_units );
		virtual void get_altitude_vector_as( double *buffer, const std::string &new_units );








		virtual void resample( double new_dz );


		virtual void as_array( double *&buffer, units_t &buffer_units );
		virtual void as_array( double *&buffer );
		virtual void as_array( NCPA::ScalarWithUnits *&buffer );

		virtual NCPA::VectorWithUnits get_vector();
		virtual NCPA::VectorWithUnits get_vector_as( units_t new_units );
		virtual NCPA::VectorWithUnits get_vector_as( const std::string &new_units );

		virtual double get( double altitude ) const;
		virtual double get_first_derivative( double altitude ) const;
		virtual double get_second_derivative( double altitude ) const;

		virtual double get_as( double altitude, units_t new_units ) const;
		virtual double get_first_derivative_as( double altitude, units_t new_units ) const;
		virtual double get_second_derivative_as( double altitude, units_t new_units ) const;

		virtual double get_as( double altitude, const std::string &new_units ) const;
		virtual double get_first_derivative_as( double altitude, const std::string &new_units ) const;
		virtual double get_second_derivative_as( double altitude, const std::string &new_units ) const;


	protected:
		NCPA::interpolator1d_t interp_type_;
		Interpolator1D *interp_;

		int check_altitude_( double z_req ) const;
		void build_splines_();
		void delete_splines_();


	};

}

#endif
