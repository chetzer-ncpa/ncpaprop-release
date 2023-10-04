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
		GSLInterpolator1D();
		GSLInterpolator1D( const gsl_interp_type *gsltype, NCPA::interpolator1d_t t );
		GSLInterpolator1D( const GSLInterpolator1D &other );
		GSLInterpolator1D( GSLInterpolator1D &&other );
		virtual ~GSLInterpolator1D();
		friend void ::swap( GSLInterpolator1D &a, GSLInterpolator1D &b );

		virtual Interpolator1D* set( size_t n, const double *x, const double *y );
		virtual Interpolator1D* set( size_t n, const double *x, const double *y_r,
				const double *y_i );
		virtual Interpolator1D* set( size_t n, const double *x,
				const std::complex<double> *y );

		virtual Interpolator1D* init();
		virtual Interpolator1D* allocate( size_t n );
		virtual Interpolator1D* ready();
		virtual void free();
		virtual bool is_ready() const;
		virtual size_t max_derivative() const;
		virtual gsl_interp_type *get_gsl_interp_type() const;
		virtual void set_gsl_interp_type( const gsl_interp_type *gsltype );

		virtual void get_interp_limits( double &xmin, double &xmax ) const;
		virtual double get_low_interp_limit() const;
		virtual double get_high_interp_limit() const;

		virtual gsl_spline *get_real_spline();
		virtual gsl_spline *get_imag_spline();
		virtual gsl_interp_accel *get_real_accel();
		virtual gsl_interp_accel *get_imag_accel();

	protected:
		bool ready_ = false;
		double minx_, maxx_;
		gsl_interp_type *interptype_ = nullptr;
		gsl_interp_accel *accel_r_ = nullptr, *accel_i_ = nullptr;
		gsl_spline *spline_r_ = nullptr, *spline_i_ = nullptr;

		void allocate_spline_( gsl_spline *&spline, gsl_interp_accel *&accel, size_t n );
		void set_splines_( size_t n, const double *x, const double *y_r,
				const double *y_i );
		void free_spline_( gsl_spline *&spline, gsl_interp_accel *&accel );

		virtual double extrapolate_f_( double x );
		virtual double extrapolate_df_( double x );
		virtual double extrapolate_d2f_( double x );
		virtual std::complex<double> extrapolate_cf_( double x );
		virtual std::complex<double> extrapolate_cdf_( double x );
		virtual std::complex<double> extrapolate_cd2f_( double x );

		virtual double interpolate_f_( double x );
		virtual double interpolate_df_( double x );
		virtual double interpolate_d2f_( double x );
		virtual std::complex<double> interpolate_cf_( double x );
		virtual std::complex<double> interpolate_cdf_( double x );
		virtual std::complex<double> interpolate_cd2f_( double x );

	};
}

#endif

#endif
