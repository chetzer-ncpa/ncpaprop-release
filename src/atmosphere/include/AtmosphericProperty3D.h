#ifndef NCPAPROP_ATMOSPHERICPROPERTY3D_H_INCLUDED
#define NCPAPROP_ATMOSPHERICPROPERTY3D_H_INCLUDED

#include <cstdlib>
#include "util.h"
#include "units.h"
#include "AtmosphericProperty1D.h"
#include "LANLInterpolation.h"

#include <vector>

typedef enum deriv_t : size_t {
	DERIVATIVE_WRT_X = 0,
	DERIVATIVE_WRT_Y = 1,
	DERIVATIVE_WRT_Z = 2 } deriv_t;

namespace NCPA {

	class AtmosphericProperty3D {

	public:
		virtual ~AtmosphericProperty3D();

		// get values
		virtual double get( double x, double y, double z ) = 0;  // one point from vector
		virtual double get( double x, double y ) = 0;            // one point from scalar
		virtual double get_derivative( double x, double y, deriv_t direction ) = 0;
		virtual double get_derivative( double x, double y, double z, deriv_t direction ) = 0;
		virtual double get_derivative( double x, double y,
			size_t order, deriv_t *directions ) = 0;
		virtual double get_derivative( double x, double y, double z,
			size_t order, deriv_t *directions ) = 0;
		virtual void as_matrix( double ***&data, size_t &nx, size_t &ny, size_t &nz ) const = 0;
		virtual void free_matrix( double ***data ) const = 0;

		// Units
		virtual units_t get_range_units() const;
		virtual units_t get_property_units() const;
		virtual units_t get_altitude_units() const;
		virtual void convert_range_units( NCPA::units_t range_units ) = 0;
		virtual void convert_property_units( NCPA::units_t new_units ) = 0;
		virtual void convert_altitude_units( NCPA::units_t new_units ) = 0;

		// limits
		virtual double x_max() const = 0;
		virtual double x_min() const = 0;
		virtual double y_max() const = 0;
		virtual double y_min() const = 0;
		virtual double z_min() const = 0;
		virtual double z_max() const = 0;
		virtual size_t x_size() const = 0;
		virtual size_t y_size() const = 0;
		virtual size_t z_size() const = 0;
		virtual void x_vector( size_t &nx, double *&x ) const = 0;
		virtual void y_vector( size_t &ny, double *&y ) const = 0;
		virtual void z_vector( size_t &nz, double *&z ) const = 0;

	protected:
		std::string key_;
		units_t r_units_;
		units_t z_units_;
		units_t f_units_;
	};

	class VectorAtmosphericProperty3D : public AtmosphericProperty3D {

	public:
		VectorAtmosphericProperty3D(
			const std::string &key,
			size_t nx, double *xvals,
			size_t ny, double *yvals,
			units_t range_units,
			NCPA::AtmosphericProperty1D ***prop_mat,
			size_t nz, double *zvals, units_t altitude_units );
		VectorAtmosphericProperty3D(
			const std::string &key, size_t nx, double *xvals, size_t ny, double *yvals,
			size_t nz, double *zvals, double ***prop_mat, NCPA::units_t range_units,
			NCPA::units_t altitude_units, NCPA::units_t property_units );
		VectorAtmosphericProperty3D( const VectorAtmosphericProperty3D &prop );
		virtual ~VectorAtmosphericProperty3D();

		// get values
		virtual double get( double x, double y, double z );  // one point from vector
		virtual double get( double x, double y );            // one point from scalar
		virtual double get_derivative( double x, double y, deriv_t direction );
		virtual double get_derivative( double x, double y, double z, deriv_t direction );
		virtual double get_derivative( double x, double y,
			size_t order, deriv_t *directions );
		virtual double get_derivative( double x, double y, double z,
			size_t order, deriv_t *directions );
		virtual void as_matrix( double ***&data, size_t &nx, size_t &ny, size_t &nz ) const;
		virtual void free_matrix( double ***data ) const;

		virtual void convert_range_units( NCPA::units_t range_units );
		virtual void convert_altitude_units( NCPA::units_t altitude_units );
		virtual void convert_property_units( NCPA::units_t new_units );

		virtual double x_max() const;
		virtual double x_min() const;
		virtual double y_max() const;
		virtual double y_min() const;
		virtual double z_max() const;
		virtual double z_min() const;
		virtual size_t x_size() const;
		virtual size_t y_size() const;
		virtual size_t z_size() const;
		virtual void x_vector( size_t &nx, double *&x ) const;
		virtual void y_vector( size_t &ny, double *&y ) const;
		virtual void z_vector( size_t &nz, double *&z ) const;

	protected:
		LANL::Spline3D spline_;

	};

	class ScalarAtmosphericProperty3D : public AtmosphericProperty3D {

	public:
		ScalarAtmosphericProperty3D(
			const std::string &key,
			size_t nx, double *xvals,
			size_t ny, double *yvals,
			units_t range_units,
			NCPA::ScalarWithUnits ***prop_mat );
		ScalarAtmosphericProperty3D(
			const std::string &key,
			size_t nx, double *xvals,
			size_t ny, double *yvals,
			double **prop_mat,
			NCPA::units_t range_units, NCPA::units_t property_units );
		ScalarAtmosphericProperty3D( const ScalarAtmosphericProperty3D &prop );
		virtual ~ScalarAtmosphericProperty3D();

		// get values
		virtual double get( double x, double y, double z );  // one point from vector
		virtual double get( double x, double y );            // one point from scalar
		virtual double get_derivative( double x, double y, deriv_t direction );
		virtual double get_derivative( double x, double y, double z, deriv_t direction );
		virtual double get_derivative( double x, double y,
			size_t order, deriv_t *directions );
		virtual double get_derivative( double x, double y, double z,
			size_t order, deriv_t *directions );
		virtual void as_matrix( double ***&data, size_t &nx, size_t &ny, size_t &nz ) const;
		virtual void free_matrix( double ***data ) const;

		virtual void convert_range_units( NCPA::units_t range_units );
		virtual void convert_property_units( NCPA::units_t new_units );
		virtual void convert_altitude_units( NCPA::units_t new_units );

		virtual double x_max() const;
		virtual double x_min() const;
		virtual double y_max() const;
		virtual double y_min() const;
		virtual double z_max() const;
		virtual double z_min() const;
		virtual size_t x_size() const;
		virtual size_t y_size() const;
		virtual size_t z_size() const;
		virtual void x_vector( size_t &nx, double *&x ) const;
		virtual void y_vector( size_t &ny, double *&y ) const;
		virtual void z_vector( size_t &nz, double *&z ) const;

	protected:
		LANL::Spline2DBicubic spline_;
	};
}

#endif