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

		virtual void set( size_t n, const double *x, const double *y );
		virtual void set( size_t n, const double *x,
				const std::complex<double> *y );
		virtual void set( size_t n, const double *x, const double *y_real,
				const double *y_imag );
		virtual void set_x( size_t n, const double *x );
		virtual void set_y( size_t n, const double *y );
		virtual void set_y( size_t n, const double *y_real, const double *y_imag );
		virtual void set_y( size_t n, const std::complex<double> *y );

		virtual void init();
		virtual void allocate( size_t n );
		virtual void ready();
		virtual void free();
		virtual bool is_ready();
		virtual size_t max_derivative() const;

		const std::string identifier() const;
		virtual double eval_f( double x );
		virtual double eval_df( double x );
		virtual double eval_ddf( double x );
		virtual double eval_dddf( double x );

		virtual std::complex<double> eval_cf( double x );
		virtual std::complex<double> eval_cdf( double x );
		virtual std::complex<double> eval_cddf( double x );
		virtual std::complex<double> eval_cdddf( double x );

//	protected:
		bool ready_ = false;
		LANL::Spline1DNatural real_spline_, imag_spline_;

		void init_spline_( LANL::Spline1DNatural &spline );
		void allocate_spline_( LANL::Spline1DNatural &spline, size_t n );
		void free_spline_( LANL::Spline1DNatural &spline );
	};
}
#endif

#endif
