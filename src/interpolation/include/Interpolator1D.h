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

		static Interpolator1D * build( interpolator1d_t t );

		Interpolator1D();
		Interpolator1D( const Interpolator1D &other );
		Interpolator1D( Interpolator1D &&other );
		virtual ~Interpolator1D() {}
		friend void ::swap( Interpolator1D &a, Interpolator1D &b );

		// interface
		virtual Interpolator1D* set( size_t n, const double *x, const double *y ) = 0;
		virtual Interpolator1D* set( size_t n, const double *x, const double *y_real,
				const double *y_imag ) = 0;
		virtual Interpolator1D* set( size_t n, const double *x,
				const std::complex<double> *y ) = 0;

		virtual double f( double x );
		virtual double df( double x );
		virtual double df( size_t n, double x );
		virtual std::complex<double> cf( double x );
		virtual std::complex<double> cdf( double x );
		virtual std::complex<double> cdf( size_t n, double x );
		virtual void is_extrapolating( bool tf );
		virtual bool is_extrapolating() const;

		virtual Interpolator1D* init() = 0;
		virtual Interpolator1D* allocate( size_t n ) = 0;
		virtual Interpolator1D* ready() = 0;
		virtual void free() = 0;
		virtual bool is_ready() = 0;
		virtual const std::string identifier() const = 0;
		virtual size_t max_derivative() const = 0;
		virtual void get_interp_limits( double &min, double &max ) const = 0;
		virtual double get_low_interp_limit() const = 0;
		virtual double get_high_interp_limit() const = 0;

	protected:
		bool extrapolate_ = true;

		virtual double eval_f_( double x );
		virtual double eval_df_( size_t n, double x );
		virtual double eval_df_( double x );
		virtual double eval_d2f_( double x );
		virtual double eval_d3f_( double x );


		virtual std::complex<double> eval_cf_( double x );
		virtual std::complex<double> eval_cdf_( size_t n, double x );
		virtual std::complex<double> eval_cdf_( double x );
		virtual std::complex<double> eval_cd2f_( double x );
		virtual std::complex<double> eval_cd3f_( double x );
	};
}


#endif
