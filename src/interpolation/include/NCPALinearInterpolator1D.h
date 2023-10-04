#ifndef NCPA__LIBINTERPOLATION_NCPALINEARINTERPOLATOR1D_H_INCLUDED
#define NCPA__LIBINTERPOLATION_NCPALINEARINTERPOLATOR1D_H_INCLUDED

#include "Interpolator1D.h"
#include <complex>

namespace NCPA { class NCPALinearInterpolator1D; }
void swap( NCPA::NCPALinearInterpolator1D &a, NCPA::NCPALinearInterpolator1D &b );

namespace NCPA {

	class NCPALinearInterpolator1D : public Interpolator1D {
	public:
		NCPALinearInterpolator1D();
		NCPALinearInterpolator1D( const NCPALinearInterpolator1D &other );
		NCPALinearInterpolator1D( NCPALinearInterpolator1D &&other );
		virtual ~NCPALinearInterpolator1D();
		friend void ::swap( NCPALinearInterpolator1D &a, NCPALinearInterpolator1D &b );
		NCPALinearInterpolator1D& operator=( NCPALinearInterpolator1D other );

		virtual Interpolator1D* set( size_t n, const double *x, const double *y );
		virtual Interpolator1D* set( size_t n, const double *x, const double *y_real,
				const double *y_imag );
		virtual Interpolator1D* set( size_t n, const double *x,
				const std::complex<double> *y );

		virtual Interpolator1D* init();
		virtual Interpolator1D* allocate( size_t n );
		virtual Interpolator1D* ready();
		virtual void free();
		virtual bool is_ready() const;
		virtual size_t max_derivative() const;

		virtual void get_interp_limits( double &xmin, double &xmax ) const;
		virtual double get_low_interp_limit() const;
		virtual double get_high_interp_limit() const;



	protected:
		size_t n_ = 0;
		double 	*x_ = nullptr,
				*yr_ = nullptr,
				*yi_ = nullptr;
		double *last_interval_ = nullptr;

		bool get_interval_( double x, size_t &i0, size_t &i1 );

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

		virtual void reset_array_( double *&x );
		virtual void allocate_array_( double *&x, size_t n );
	};
}


#endif
