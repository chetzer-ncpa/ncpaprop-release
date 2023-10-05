#include "Interpolator1D.h"

#ifndef NCPA__LIBINTERPOLATION_LANLNATURALCUBICSPLINEINTERPOLATOR1D_H_INCLUDED
#define NCPA__LIBINTERPOLATION_LANLNATURALCUBICSPLINEINTERPOLATOR1D_H_INCLUDED

#ifdef HAVE_LANL_INTERPOLATION_LIBRARY
#include "LANLInterpolation.h"
#include <complex>

namespace NCPA { class LANLNaturalCubicSplineInterpolator1D; }
void swap( NCPA::LANLNaturalCubicSplineInterpolator1D &a, NCPA::LANLNaturalCubicSplineInterpolator1D &b );

namespace NCPA {

	class LANLNaturalCubicSplineInterpolator1D : public Interpolator1D {
	public:
		LANLNaturalCubicSplineInterpolator1D();
		LANLNaturalCubicSplineInterpolator1D( const LANLNaturalCubicSplineInterpolator1D &other );
		LANLNaturalCubicSplineInterpolator1D( LANLNaturalCubicSplineInterpolator1D &&other );
		virtual ~LANLNaturalCubicSplineInterpolator1D();
		friend void ::swap( LANLNaturalCubicSplineInterpolator1D &a, LANLNaturalCubicSplineInterpolator1D &b );
		virtual Interpolator1D* clone() const;

		virtual Interpolator1D* set( size_t n, const double *x, const double *y );
		virtual Interpolator1D* set( size_t n, const double *x,
				const std::complex<double> *y );
		virtual Interpolator1D* set( size_t n, const double *x, const double *y_real,
				const double *y_imag );

		virtual void set_x( size_t n, const double *x );
		virtual void set_y( size_t n, const double *y );
		virtual void set_y( size_t n, const double *y_real, const double *y_imag );
		virtual void set_y( size_t n, const std::complex<double> *y );

		virtual Interpolator1D* init();
		virtual Interpolator1D* allocate( size_t n );
		virtual Interpolator1D* ready();
		virtual void free();
		virtual bool is_ready() const;
		virtual size_t max_derivative() const;

		virtual void get_interp_limits( double &xmin, double &xmax ) const;
		virtual double get_low_interp_limit() const;
		virtual double get_high_interp_limit() const;

		virtual LANL::Spline1DNatural *get_real_spline();
		virtual LANL::Spline1DNatural *get_imag_spline();

	protected:
		bool ready_ = false;
		double minx_, maxx_;
		LANL::Spline1DNatural real_spline_, imag_spline_;

		void init_spline_( LANL::Spline1DNatural &spline );
		void allocate_spline_( LANL::Spline1DNatural &spline, size_t n );
		void free_spline_( LANL::Spline1DNatural &spline );

		virtual double interpolate_f_( double x );
		virtual double interpolate_df_( double x );
		virtual double interpolate_d2f_( double x );
		virtual double interpolate_d3f_( double x );
		virtual std::complex<double> interpolate_cf_( double x );
		virtual std::complex<double> interpolate_cdf_( double x );
		virtual std::complex<double> interpolate_cd2f_( double x );
		virtual std::complex<double> interpolate_cd3f_( double x );

		virtual double extrapolate_f_( double x );
		virtual double extrapolate_df_( double x );
		virtual double extrapolate_d2f_( double x );
		virtual double extrapolate_d3f_( double x );
		virtual std::complex<double> extrapolate_cf_( double x );
		virtual std::complex<double> extrapolate_cdf_( double x );
		virtual std::complex<double> extrapolate_cd2f_( double x );
		virtual std::complex<double> extrapolate_cd3f_( double x );
	};
}
#endif

#endif
