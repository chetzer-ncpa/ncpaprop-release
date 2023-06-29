#ifndef NCPA__LIBINTERPOLATION_GSLINTERPOLATOR1D_H_INCLUDED
#define NCPA__LIBINTERPOLATION_GSLINTERPOLATOR1D_H_INCLUDED

#include "Interpolator1D.h"

#ifdef HAVE_GSL_INTERPOLATION_LIBRARY

#include "gsl/gsl_spline.h"
#include <complex>

namespace NCPA { class GSLInterpolator1D; }
void swap( NCPA::GSLInterpolator1D &a, NCPA::GSLInterpolator1D &b );

namespace NCPA {

	class GSLInterpolator1D : public Interpolator1D {
	public:
		GSLInterpolator1D( const gsl_interp_type *gsltype );
		GSLInterpolator1D( const GSLInterpolator1D &other );
		GSLInterpolator1D( GSLInterpolator1D &&other );
		virtual ~GSLInterpolator1D();
		friend void ::swap( GSLInterpolator1D &a, GSLInterpolator1D &b );

		virtual void set( size_t n, const double *x, const double *y );
		virtual void set( size_t n, const double *x, const double *y_r,
				const double *y_i );
		virtual void set( size_t n, const double *x,
				const std::complex<double> *y );

		virtual void init();
		virtual void allocate( size_t n );
		virtual void ready();
		virtual void free();
		virtual bool is_ready();
		virtual size_t max_derivative() const;
		virtual gsl_interp_type *get_gsl_interp_type() const;

		const std::string identifier() const;
		virtual double eval_f( double x );
		virtual double eval_df( double x );
		virtual double eval_ddf( double x );

		virtual std::complex<double> eval_cf( double x );
		virtual std::complex<double> eval_cdf( double x );
		virtual std::complex<double> eval_cddf( double x );

//	protected:
		bool ready_ = false;
		gsl_interp_type *interptype_;
		gsl_interp_accel *accel_r_ = nullptr, *accel_i_ = nullptr;
		gsl_spline *spline_r_ = nullptr, *spline_i_ = nullptr;

		void allocate_spline_( gsl_spline *&spline, gsl_interp_accel *&accel, size_t n );
		void set_splines_( size_t n, const double *x, const double *y_r,
				const double *y_i );
		void free_spline_( gsl_spline *&spline, gsl_interp_accel *&accel );
	};
}

#endif

#endif
