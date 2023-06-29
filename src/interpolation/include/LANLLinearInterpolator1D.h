#include "Interpolator1D.h"

#ifndef NCPA__LIBINTERPOLATION_LANLLINEARINTERPOLATOR1D_H_INCLUDED
#define NCPA__LIBINTERPOLATION_LANLLINEARINTERPOLATOR1D_H_INCLUDED

#ifdef HAVE_LANL_INTERPOLATION_LIBRARY

#include "LANLInterpolation.h"
#include <complex>

namespace NCPA { class LANLLinearInterpolator1D; }
void swap( NCPA::LANLLinearInterpolator1D &a, NCPA::LANLLinearInterpolator1D &b );

namespace NCPA {

	class LANLLinearInterpolator1D : public Interpolator1D {
	public:
		LANLLinearInterpolator1D();
		LANLLinearInterpolator1D( const LANLLinearInterpolator1D &other );
		LANLLinearInterpolator1D( LANLLinearInterpolator1D &&other );
		virtual ~LANLLinearInterpolator1D();
		friend void ::swap( LANLLinearInterpolator1D &a, LANLLinearInterpolator1D &b );

		virtual void set( size_t n, const double *x, const double *y );
		virtual void set( size_t n, const double *x, const double *y_real,
				const double *y_imag );
		virtual void set( size_t n, const double *x,
				const std::complex<double> *y );
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

		virtual std::complex<double> eval_cf( double x );
		virtual std::complex<double> eval_cdf( double x );

//	protected:
		bool ready_ = false;
		LANL::Spline1DLinear real_spline_, imag_spline_;

		void init_spline_( LANL::Spline1DLinear &spline );
		void allocate_spline_( LANL::Spline1DLinear &spline, size_t n );
		void free_spline_( LANL::Spline1DLinear &spline );
	};
}

#endif

#endif
