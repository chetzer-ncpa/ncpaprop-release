#ifndef NCPA__INTERPOLATION_GSLPERIODICINTERPOLATOR_H_INCLUDED
#define NCPA__INTERPOLATION_GSLPERIODICINTERPOLATOR_H_INCLUDED

#include "Interpolator1D.h"

#ifdef HAVE_GSL_INTERPOLATION_LIBRARY

#include "gsl/gsl_spline.h"
#include "Interpolator1D.h"
#include "GSLInterpolator1D.h"
#include <complex>

namespace NCPA { class GSLPeriodicInterpolator1D; }
void swap( NCPA::GSLPeriodicInterpolator1D &a, NCPA::GSLPeriodicInterpolator1D &b );

namespace NCPA {

	class GSLPeriodicInterpolator1D : public GSLInterpolator1D {
	public:
		GSLPeriodicInterpolator1D( const gsl_interp_type *gsltype,
				NCPA::interpolator1d_t ncpatype );
		GSLPeriodicInterpolator1D( const GSLPeriodicInterpolator1D &other );
		GSLPeriodicInterpolator1D( GSLPeriodicInterpolator1D &&other );
		virtual ~GSLPeriodicInterpolator1D();
		friend void ::swap( GSLPeriodicInterpolator1D &a, GSLPeriodicInterpolator1D &b );

		virtual double f( double x );
		virtual double df( double x );
		virtual double df( size_t n, double x );
		virtual std::complex<double> cf( double x );
		virtual std::complex<double> cdf( double x );
		virtual std::complex<double> cdf( size_t n, double x );
	};

}
#endif
#endif
