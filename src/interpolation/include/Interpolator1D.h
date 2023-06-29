#ifndef NCPA__INTERPOLATOR1D_H_INCLUDED
#define NCPA__INTERPOLATOR1D_H_INCLUDED

#include <complex>
#include <vector>

#if __has_include("LANLInterpolation.h")
#define HAVE_LANL_INTERPOLATION_LIBRARY
#endif

#if __has_include("gsl/gsl_spline.h")
#define HAVE_GSL_INTERPOLATION_LIBRARY
#endif

namespace NCPA { class Interpolator1D; }
void swap( NCPA::Interpolator1D &a, NCPA::Interpolator1D &b );

namespace NCPA {

	enum class interpolator1d_t : unsigned int {
		NEAREST_NEIGHBOR = 0,
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
		Interpolator1D();
		Interpolator1D( const Interpolator1D &other );
		Interpolator1D( Interpolator1D &&other );
		virtual ~Interpolator1D() {}
		friend void ::swap( Interpolator1D &a, Interpolator1D &b );

		// interface
		virtual void set( size_t n, const double *x, const double *y ) = 0;
		virtual void set( size_t n, const double *x, const double *y_real,
				const double *y_imag ) = 0;
		virtual void set( size_t n, const double *x,
				const std::complex<double> *y ) = 0;

		virtual void init() = 0;
		virtual void allocate( size_t n ) = 0;
		virtual void ready() = 0;
		virtual void free() = 0;

		virtual bool is_ready() = 0;
		virtual const std::string identifier() const = 0;
		virtual size_t max_derivative() const = 0;

		virtual double eval_f( double x );
		virtual double eval_df( size_t n, double x );
		virtual double eval_df( double x );
		virtual double eval_ddf( double x );
		virtual double eval_dddf( double x );


		virtual std::complex<double> eval_cf( double x );
		virtual std::complex<double> eval_cdf( size_t n, double x );
		virtual std::complex<double> eval_cdf( double x );
		virtual std::complex<double> eval_cddf( double x );
		virtual std::complex<double> eval_cdddf( double x );

		static Interpolator1D * build( interpolator1d_t t );
	};
}


#endif
