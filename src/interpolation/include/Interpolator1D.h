#ifndef NCPA__INTERPOLATOR1D_H_INCLUDED
#define NCPA__INTERPOLATOR1D_H_INCLUDED

#include <complex>
#include <vector>

#if __has_include("LANLInterpolation.h")
#define HAVE_LANL_INTERPOLATION_LIBRARY
#include "LANLInterpolation.h"
#endif

#if __has_include("gsl/gsl_spline.h")
#define HAVE_GSL_INTERPOLATION_LIBRARY
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_version.h"
#endif

namespace NCPA { class Interpolator1D; }
void swap( NCPA::Interpolator1D &a, NCPA::Interpolator1D &b );

namespace NCPA {

	enum class interpolator1d_t : unsigned int {
		INVALID = 0,
		NCPA_1D_NEAREST_NEIGHBOR,
		NCPA_1D_LINEAR,
		LANL_1D_LINEAR,
		LANL_1D_NATURAL_CUBIC,
		GSL_1D_LINEAR,
		GSL_1D_POLYNOMIAL,
		GSL_1D_CUBIC,
		GSL_1D_CUBIC_PERIODIC,
		GSL_1D_AKIMA,
		GSL_1D_AKIMA_PERIODIC,
		GSL_1D_STEFFEN,
	};

	class Interpolator1D {
	public:

		static Interpolator1D * build( interpolator1d_t t );
		static bool can_build( interpolator1d_t t );
		static std::string as_string( interpolator1d_t t );

		Interpolator1D();
		Interpolator1D( interpolator1d_t t );
		Interpolator1D( const Interpolator1D &other );
		Interpolator1D( Interpolator1D &&other );
		virtual ~Interpolator1D() {}
		friend void ::swap( Interpolator1D &a, Interpolator1D &b );
		virtual Interpolator1D* clone() const = 0;

		// interface
		virtual Interpolator1D* set( size_t n, const double *x, const double *y ) = 0;
		virtual Interpolator1D* set( size_t n, const double *x, const double *y_real,
				const double *y_imag ) = 0;
		virtual Interpolator1D* set( size_t n, const double *x,
				const std::complex<double> *y ) = 0;

		virtual Interpolator1D* set( std::vector<double> x, std::vector<double> y );
		virtual Interpolator1D* set( std::vector<double> x, std::vector<std::complex<double>> y );

		virtual double f( double x );
		virtual double df( double x );
		virtual double df( size_t n, double x );
		virtual std::complex<double> cf( double x );
		virtual std::complex<double> cdf( double x );
		virtual std::complex<double> cdf( size_t n, double x );
		virtual void is_extrapolating( bool tf );
		virtual bool is_extrapolating() const;
		virtual const std::string identifier() const;
		virtual const interpolator1d_t type() const;

		virtual Interpolator1D* init() = 0;
		virtual Interpolator1D* allocate( size_t n ) = 0;
		virtual Interpolator1D* ready() = 0;
		virtual void free() = 0;
		virtual bool is_ready() const = 0;
		virtual size_t max_derivative() const = 0;
		virtual void get_interp_limits( double &min, double &max ) const = 0;
		virtual double get_low_interp_limit() const = 0;
		virtual double get_high_interp_limit() const = 0;

	protected:
		bool extrapolating_ = true;
		interpolator1d_t type_;

		// convenience functions for the default linear extrapolation,
		// i.e. for f, linearly extrapolate based on the first derivatives
		// at the end points; for df, report the first derivative at the
		// appropriate endpoint, for higher derivatives return 0.0.
		virtual double linear_extrapolate_f_( double x ) final;
		virtual double linear_extrapolate_df_( double x ) final;
		virtual double linear_extrapolate_d2f_( double x ) final;
		virtual double linear_extrapolate_d3f_( double x ) final;
		virtual std::complex<double> linear_extrapolate_cf_( double x ) final;
		virtual std::complex<double> linear_extrapolate_cdf_( double x ) final;
		virtual std::complex<double> linear_extrapolate_cd2f_( double x ) final;
		virtual std::complex<double> linear_extrapolate_cd3f_( double x ) final;


		// These functions are to be overridden by subclasses.  If not
		// overridden, they will throw domain_error because not applicable
		// or not implemented for that type of interpolator
		virtual double extrapolate_( size_t n, double x );
		virtual double extrapolate_f_( double x );
		virtual double extrapolate_df_( double x );
		virtual double extrapolate_d2f_( double x );
		virtual double extrapolate_d3f_( double x );

		virtual std::complex<double> extrapolate_c_( size_t n, double x );
		virtual std::complex<double> extrapolate_cf_( double x );
		virtual std::complex<double> extrapolate_cdf_( double x );
		virtual std::complex<double> extrapolate_cd2f_( double x );
		virtual std::complex<double> extrapolate_cd3f_( double x );

		virtual double interpolate_( size_t n, double x );
		virtual double interpolate_f_( double x );
		virtual double interpolate_df_( double x );
		virtual double interpolate_d2f_( double x );
		virtual double interpolate_d3f_( double x );

		virtual std::complex<double> interpolate_c_( size_t n, double x );
		virtual std::complex<double> interpolate_cf_( double x );
		virtual std::complex<double> interpolate_cdf_( double x );
		virtual std::complex<double> interpolate_cd2f_( double x );
		virtual std::complex<double> interpolate_cd3f_( double x );
	};
}


#endif
