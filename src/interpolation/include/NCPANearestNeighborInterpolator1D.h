#ifndef NCPA__INTERPOLATION_NCPANEARESTNEIGHBORINTERPOLATOR1D_H_INCLUDED
#define NCPA__INTERPOLATION_NCPANEARESTNEIGHBORINTERPOLATOR1D_H_INCLUDED

#include "Interpolator1D.h"
#include <complex>

namespace NCPA { class NCPANearestNeighborInterpolator1D; }
void swap( NCPA::NCPANearestNeighborInterpolator1D &a, NCPA::NCPANearestNeighborInterpolator1D &b );

namespace NCPA {

	class NCPANearestNeighborInterpolator1D : public Interpolator1D {
	public:
		NCPANearestNeighborInterpolator1D();
		NCPANearestNeighborInterpolator1D( const NCPANearestNeighborInterpolator1D &other );
		NCPANearestNeighborInterpolator1D( NCPANearestNeighborInterpolator1D &&other );
		virtual ~NCPANearestNeighborInterpolator1D();
		friend void ::swap( NCPANearestNeighborInterpolator1D &a, NCPANearestNeighborInterpolator1D &b );
		NCPANearestNeighborInterpolator1D& operator=( NCPANearestNeighborInterpolator1D other );

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
		virtual void get_interp_limits( double &min, double &max ) const;
		virtual double get_low_interp_limit() const;
		virtual double get_high_interp_limit() const;
//		virtual bool is_extrapolating() const;

	protected:
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

		size_t n_ = 0;
		double *x_ = nullptr, *yr_ = nullptr, *yi_ = nullptr;
	};

}
#endif
