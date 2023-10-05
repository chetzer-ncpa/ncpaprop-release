#include "NCPAInterpolation.h"
#include "NCPACommon.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <complex>
#include <iostream>

using namespace std;
using namespace NCPA;
using namespace testing;

#ifndef EXPECT_DOUBLE_ARRAY_NEAR
#define EXPECT_DOUBLE_ARRAY_NEAR(N,A,Ex,T) for (int i=0; i<N; i++) { EXPECT_NEAR(A[i],Ex[i],T); }
#endif

#ifndef EXPECT_DOUBLE_ARRAY_EQ
#define EXPECT_DOUBLE_ARRAY_EQ(N,A,Ex) for (int i=0; i<N; i++) { EXPECT_DOUBLE_EQ(A[i],Ex[i]); }
#endif


// test fixture
class GSLInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		double phase = M_PI / 4.0;
		for (size_t i = 0; i < 9; i++) {
			// test function for nonlinear aperiodic interpolators is
			// a cubic polynomial
			// y = 2x**3 + 1.25x**2 - 1.5x + 2
			// dy/dx = 6x**2 + 2.5x - 1.5
			// d2y/dx2 = 12x + 2.5
			n_rvals[i] = evalpoly( 4, coeffs, n_xvals[i] );
			n_ivals[i] = 0.5 * n_rvals[i];
			if (i > 0) {
				n_interp_xvals[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_interp_rvals[i-1] = evalpoly( 4, coeffs, n_interp_xvals[i-1] );
				n_interp_ivals[i-1] = 0.5 * n_interp_rvals[i-1];
				n_deriv1_xvals[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_deriv1_rvals[i-1] = evalpoly( 3, coeffs1d, n_deriv1_xvals[i-1] );
				n_deriv1_ivals[i-1] = 0.5 * n_deriv1_rvals[i-1];
				n_deriv2_xvals[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_deriv2_rvals[i-1] = evalpoly( 2, coeffs2d, n_deriv2_xvals[i-1] );
				n_deriv2_ivals[i-1] = 0.5 * n_deriv2_rvals[i-1];
			}

		}
		linear = Interpolator1D::build(interpolator1d_t::GSL_1D_LINEAR);
		polynomial = Interpolator1D::build(interpolator1d_t::GSL_1D_POLYNOMIAL);
		cubic = Interpolator1D::build(interpolator1d_t::GSL_1D_CUBIC);
		akima = Interpolator1D::build(interpolator1d_t::GSL_1D_AKIMA);
		steffen = Interpolator1D::build(interpolator1d_t::GSL_1D_STEFFEN);

//		cout << "Setup: Linear" << endl;
//		for (size_t i = 0; i < 5; i++) {
//			cout << "xvals[i]=" << l_xvals[i] << endl;
//			cout << "yvals_r[i]=" << l_rvals[i] << endl;
//			cout << "yvals_i[i]=" << l_ivals[i] << endl;
//		}
//		cout << "Linear setup done" << endl;
		linear->set( 5, l_xvals, l_rvals, l_ivals )->ready();

//		cout << "Setup: Nonlinear" << endl;
//		for (size_t i = 0; i < 9; i++) {
//			cout << "xvals[i]=" << n_xvals[i] << endl;
//			cout << "yvals_r[i]=" << n_rvals[i] << endl;
//			cout << "yvals_i[i]=" << n_ivals[i] << endl;
//		}
//		cout << "Nonlinear setup done" << endl;

		polynomial->set( 9, n_xvals, n_rvals, n_ivals )->ready();
		cubic->set( 9, n_xvals, n_rvals, n_ivals )->ready();
		akima->set( 9, n_xvals, n_rvals, n_ivals )->ready();
		steffen->set( 9, n_xvals, n_rvals, n_ivals )->ready();

		double2complex( 5, l_rvals, l_ivals, l_cvals );
	}

	// If there's any cleanup other than normal destructors
//	void TearDown() override {
//		delete interp1;
//		delete c_interp1;
//	}

	// class members here
	Interpolator1D *linear, *polynomial, *cubic,
					*akima, *steffen;
	GSLInterpolator1D *gPtr1, *gPtr2;
	double l_tolerance = 1e-3;
	double n_tolerance = 0.15;

	// vectors for testing linear interpolators
	double l_xvals[5] = { 1, 2, 3, 4, 5 };
	double l_rvals[5] = { 2, 4, 8, 16, 32 };
	double l_ivals[5] = { 2, 1, 0, -2, -4 };
	double l_interp_x[4] = { 1.5, 2.5, 3.5, 4.5 };
	double l_interp_r[4] = { 3, 6, 12, 24 };
	double l_interp_i[4] = { 1.5, 0.5, -1, -3 };
	double l_derivs_x[4] = { 1.5, 2.5, 3.5, 4.5 };
	double l_derivs_r[4] = { 2, 4, 8, 16 };
	double l_derivs_i[4] = { -1, -1, -2, -2 };
	double l_derivs2[4] = { 0, 0, 0, 0 };
	complex<double> l_cvals[5];

	// vectors for testing nonlinear interpolators
	double n_xvals[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double coeffs[4] = { 2.0, -1.5, 1.25, 2.0 };
	double coeffs1d[3] = { -1.5, 2.5, 6.0 };
	double coeffs2d[2] = { 2.5, 12.0 };

	// to calculate
	double n_rvals[9], n_ivals[9];
	complex<double> n_cvals[9];
	double n_interp_xvals[8], n_interp_rvals[8], n_interp_ivals[8];
	double n_deriv1_xvals[8], n_deriv1_rvals[8], n_deriv1_ivals[8];
	double n_deriv2_xvals[8], n_deriv2_rvals[8], n_deriv2_ivals[8];

	double zvals5[5] = { 0, 0, 0, 0, 0 };
	double zvals9[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
};


// Tests that use the fixture
TEST_F(GSLInterpolator1DTest,CopyConstructorCreatesCopy) {
	complex<double> testval = linear->cf( 2.1 );
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	GSLInterpolator1D linear2( *gPtr1 );

	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_real_spline()->x, gPtr1->get_real_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_real_spline()->y, gPtr1->get_real_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_imag_spline()->x, gPtr1->get_imag_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_imag_spline()->y, gPtr1->get_imag_spline()->y );
	EXPECT_EQ( linear2.get_real_accel()->cache, gPtr1->get_real_accel()->cache );
	EXPECT_EQ( linear2.get_real_accel()->miss_count, gPtr1->get_real_accel()->miss_count );
	EXPECT_EQ( linear2.get_real_accel()->hit_count, gPtr1->get_real_accel()->hit_count );
	EXPECT_EQ( linear2.get_imag_accel()->cache, gPtr1->get_imag_accel()->cache );
	EXPECT_EQ( linear2.get_imag_accel()->miss_count, gPtr1->get_imag_accel()->miss_count );
	EXPECT_EQ( linear2.get_imag_accel()->hit_count, gPtr1->get_imag_accel()->hit_count );
	EXPECT_TRUE( linear2.is_ready() );
}

TEST_F(GSLInterpolator1DTest,SwapWorksProperly) {
	double testval = steffen->f( 2.5 );
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	gPtr2 = static_cast<GSLInterpolator1D *>( steffen );
	GSLInterpolator1D linear2( *gPtr1 ),
					  steffen2( *gPtr2 );
	EXPECT_EQ( linear2.get_real_accel()->cache, gPtr1->get_real_accel()->cache );
	EXPECT_EQ( steffen2.get_real_accel()->cache, gPtr2->get_real_accel()->cache );
	EXPECT_EQ( linear2.get_real_accel()->miss_count, gPtr1->get_real_accel()->miss_count );
	EXPECT_EQ( steffen2.get_real_accel()->miss_count, gPtr2->get_real_accel()->miss_count );
	EXPECT_EQ( linear2.get_real_accel()->hit_count, gPtr1->get_real_accel()->hit_count );
	EXPECT_EQ( steffen2.get_real_accel()->hit_count, gPtr2->get_real_accel()->hit_count );

	swap( linear2, steffen2 );

	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_real_spline()->x, gPtr2->get_real_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_real_spline()->y, gPtr2->get_real_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_imag_spline()->x, gPtr2->get_imag_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.get_imag_spline()->y, gPtr2->get_imag_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.get_real_spline()->x, gPtr1->get_real_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.get_real_spline()->y, gPtr1->get_real_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.get_imag_spline()->x, gPtr1->get_imag_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.get_imag_spline()->y, gPtr1->get_imag_spline()->y );
	EXPECT_EQ( linear2.identifier(), gPtr2->identifier() );
	EXPECT_EQ( steffen2.identifier(), gPtr1->identifier() );
	EXPECT_EQ( linear2.get_real_accel()->cache, gPtr2->get_real_accel()->cache );
	EXPECT_EQ( steffen2.get_real_accel()->cache, gPtr1->get_real_accel()->cache );
	EXPECT_EQ( linear2.get_real_accel()->miss_count, gPtr2->get_real_accel()->miss_count );
	EXPECT_EQ( steffen2.get_real_accel()->miss_count, gPtr1->get_real_accel()->miss_count );
	EXPECT_EQ( linear2.get_real_accel()->hit_count, gPtr2->get_real_accel()->hit_count );
	EXPECT_EQ( steffen2.get_real_accel()->hit_count, gPtr1->get_real_accel()->hit_count );
}

TEST_F(GSLInterpolator1DTest,CloneWorksProperly) {
	Interpolator1D *clone = linear->clone();
	linear->free();
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( clone->f(l_xvals[i]), l_rvals[i] );
	}
	delete clone;
	clone = polynomial->clone();
	polynomial->free();
	for (size_t i = 0; i < 9; i++) {
		EXPECT_DOUBLE_EQ( clone->f(n_xvals[i]), n_rvals[i] );
	}
	delete clone;
	clone = steffen->clone();
	steffen->free();
	for (size_t i = 0; i < 9; i++) {
		EXPECT_DOUBLE_EQ( clone->f(n_xvals[i]), n_rvals[i] );
	}
	delete clone;
	clone = akima->clone();
	akima->free();
	for (size_t i = 0; i < 9; i++) {
		EXPECT_DOUBLE_EQ( clone->f(n_xvals[i]), n_rvals[i] );
	}
	delete clone;
	clone = cubic->clone();
	cubic->free();
	for (size_t i = 0; i < 9; i++) {
		EXPECT_DOUBLE_EQ( clone->f(n_xvals[i]), n_rvals[i] );
	}
	delete clone;
}

TEST_F(GSLInterpolator1DTest,RealSetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->y, l_rvals );
	gPtr1->set( 9, n_xvals, n_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->y, n_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->y, zvals9 );

}

TEST_F(GSLInterpolator1DTest,RealAndImaginarySetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->y, l_rvals );
	gPtr1->set( 9, n_xvals, n_rvals, n_ivals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->y, n_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->y, n_ivals );
}

TEST_F(GSLInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLInterpolator1D *>( steffen );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->y, n_rvals );
	gPtr1->set( 5, l_xvals, l_cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->y, l_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_imag_spline()->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_imag_spline()->y, l_ivals );
}

TEST_F(GSLInterpolator1DTest,InitSetsPointersToNull) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	gPtr1->init();
	EXPECT_EQ( gPtr1->get_real_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_real_accel(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_accel(), nullptr );
}

TEST_F(GSLInterpolator1DTest,AllocateSetsPointersToValid) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	gPtr1->init();
	gPtr1->allocate(15);
	EXPECT_NE( gPtr1->get_real_spline(), nullptr );
	EXPECT_NE( gPtr1->get_imag_spline(), nullptr );
	EXPECT_NE( gPtr1->get_real_accel(), nullptr );
	EXPECT_NE( gPtr1->get_imag_accel(), nullptr );
}

TEST_F(GSLInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	gPtr1->init();
}

TEST_F(GSLInterpolator1DTest,FreeSetsPointersToNull) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	cubic->free();
	EXPECT_EQ( gPtr1->get_real_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_real_accel(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_accel(), nullptr );
}

TEST_F(GSLInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( linear->is_ready() );
}

TEST_F(GSLInterpolator1DTest,GetInterpLimitsReturnsCorrectLimits) {
	double x1, x2;
	steffen->get_interp_limits(x1, x2);
	EXPECT_DOUBLE_EQ( x1, 1.0 );
	EXPECT_DOUBLE_EQ( x2, 9.0 );
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectRealValue) {
	EXPECT_DOUBLE_EQ( linear->f(-1.0), -2.0 );

	// nonlinear aperiodic should all work the same:
	// f(1) = 3.75
	double x0 = 1.0;
	double f0 = polynomial->f(x0);
	double stepsize = -1.0;
	size_t steps = 1;
	double x = x0 + (double)steps * stepsize;
	EXPECT_DOUBLE_EQ( polynomial->f(x), f0 - (double)steps * polynomial->df(x0) );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( cubic->f(x), f0 - (double)steps * cubic->df(x0) );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( akima->f(x), f0 - (double)steps * akima->df(x0) );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( steffen->f(x), f0 - (double)steps * steffen->df(x0) );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectRealFirstDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	EXPECT_DOUBLE_EQ( linear->df(-1.0), linear->df(1.0) );
	EXPECT_DOUBLE_EQ( polynomial->df(0.0), polynomial->df(1.0) );
	EXPECT_DOUBLE_EQ( cubic->df(-1.0), cubic->df(1.0) );
	EXPECT_DOUBLE_EQ( akima->df(-2.0), akima->df(1.0) );
	EXPECT_DOUBLE_EQ( steffen->df(-3.0), steffen->df(1.0) );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectRealSecondDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	EXPECT_DOUBLE_EQ( linear->df(2,-1.0), 0.0 );
	EXPECT_DOUBLE_EQ( polynomial->df(2,0.0), 0.0 );
	EXPECT_DOUBLE_EQ( cubic->df(2,-1.0), 0.0 );
	EXPECT_DOUBLE_EQ( akima->df(2,-2.0), 0.0 );
	EXPECT_DOUBLE_EQ( steffen->df(2,-3.0), 0.0 );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectRealValue) {
	EXPECT_DOUBLE_EQ( linear->f(7.0), 64.0 );

	// nonlinear aperiodic should all work the same:
	// f(9) = 1547.75
	double x0 = 9.0;
	double f0 = polynomial->f(x0);
	double stepsize = 1.0;
	size_t steps = 1;
	double x = x0 + (double)steps * stepsize;
	EXPECT_DOUBLE_EQ( polynomial->f(x), f0 + (double)steps * polynomial->df(x0) );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( cubic->f(x), f0 + (double)steps * cubic->df(x0) );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( akima->f(x), f0 + (double)steps * akima->df(x0) );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( steffen->f(x), f0 + (double)steps * steffen->df(x0) );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectRealFirstDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	double x0 = linear->get_high_interp_limit();
	double xtest = 7.0;
	EXPECT_DOUBLE_EQ( linear->df(xtest), linear->df(x0) );

	x0 = polynomial->get_high_interp_limit();
	xtest = x0 + 2.0;
	EXPECT_DOUBLE_EQ( polynomial->df(xtest), polynomial->df(x0) );
	EXPECT_DOUBLE_EQ( cubic->df(xtest+1.0), cubic->df(x0) );
	EXPECT_DOUBLE_EQ( akima->df(xtest+1.5), akima->df(x0) );
	EXPECT_DOUBLE_EQ( steffen->df(xtest+2.5), steffen->df(x0) );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectRealSecondDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	double x0 = linear->get_high_interp_limit();
	double xtest = 7.0;
	EXPECT_DOUBLE_EQ( linear->df(2,-1.0), 0.0 );

	x0 = polynomial->get_high_interp_limit();
	xtest = x0 + 2.0;
	EXPECT_DOUBLE_EQ( polynomial->df(2,xtest), 0.0 );
	EXPECT_DOUBLE_EQ( cubic->df(2,xtest+1.0), 0.0 );
	EXPECT_DOUBLE_EQ( akima->df(2,xtest+3.0), 0.0 );
	EXPECT_DOUBLE_EQ( steffen->df(2,xtest+5.24), 0.0 );

	// test periodic in a minute
}


TEST_F(GSLInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectComplexValue) {
	EXPECT_DOUBLE_EQ( linear->cf(-1.0).real(), -2.0 );
	EXPECT_DOUBLE_EQ( linear->cf(-1.0).imag(), 4.0 );

	// nonlinear aperiodic should all work the same:
	// f(1) = 3.75
	double x0 = 1.0;
	complex<double> f0 = polynomial->cf(x0);
	double stepsize = -1.0;
	size_t steps = 1;

	double x = x0 + (double)steps * stepsize;
	EXPECT_DOUBLE_EQ( polynomial->cf(x).real(), (f0 - (double)steps * polynomial->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( polynomial->cf(x).imag(), (f0 - (double)steps * polynomial->cdf(x0)).imag() );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( cubic->cf(x).real(), (f0 - (double)steps * cubic->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( cubic->cf(x).imag(), (f0 - (double)steps * cubic->cdf(x0)).imag() );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( akima->cf(x).real(), (f0 - (double)steps * akima->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( akima->cf(x).imag(), (f0 - (double)steps * akima->cdf(x0)).imag() );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( steffen->cf(x).real(), (f0 - (double)steps * steffen->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( steffen->cf(x).imag(), (f0 - (double)steps * steffen->cdf(x0)).imag() );
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectComplexFirstDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	double x0 = linear->get_low_interp_limit();
	double step = -2.0;
	EXPECT_DOUBLE_EQ( linear->cdf(x0 + step).real(), linear->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( linear->cdf(x0 + step).imag(), linear->cdf(x0).imag() );

	x0 = polynomial->get_low_interp_limit();
	EXPECT_DOUBLE_EQ( polynomial->cdf(x0 + step).real(), polynomial->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( cubic->cdf(x0 + step-1).real(), cubic->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( akima->cdf(x0 + step-2).real(), akima->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( steffen->cdf(x0 + step-3).real(), steffen->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( polynomial->cdf(x0 + step).imag(), polynomial->cdf(x0).imag() );
	EXPECT_DOUBLE_EQ( cubic->cdf(x0 + step-1).imag(), cubic->cdf(x0).imag() );
	EXPECT_DOUBLE_EQ( akima->cdf(x0 + step-2).imag(), akima->cdf(x0).imag() );
	EXPECT_DOUBLE_EQ( steffen->cdf(x0 + step-3).imag(), steffen->cdf(x0).imag() );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectComplexSecondDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	double xtest = linear->get_low_interp_limit() - 5.0;
	EXPECT_DOUBLE_EQ( linear->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( polynomial->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( cubic->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( akima->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( steffen->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( linear->cdf(2,xtest).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( polynomial->cdf(2,xtest).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( cubic->cdf(2,xtest).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( akima->cdf(2,xtest).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( steffen->cdf(2,xtest).imag(), 0.0 );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectComplexValue) {
	EXPECT_DOUBLE_EQ( linear->cf(7.0).real(), 64.0 );
	EXPECT_DOUBLE_EQ( linear->cf(7.0).imag(), -8.0 );

	// nonlinear aperiodic should all work the same:
	// f(1) = 3.75
	double x0 = polynomial->get_high_interp_limit();
	complex<double> f0 = polynomial->cf(x0);
	double stepsize = 1.0;
	size_t steps = 1;

	double x = x0 + (double)steps * stepsize;
	EXPECT_DOUBLE_EQ( polynomial->cf(x).real(), (f0 + (double)steps * polynomial->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( polynomial->cf(x).imag(), (f0 + (double)steps * polynomial->cdf(x0)).imag() );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( cubic->cf(x).real(), (f0 + (double)steps * cubic->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( cubic->cf(x).imag(), (f0 + (double)steps * cubic->cdf(x0)).imag() );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( akima->cf(x).real(), (f0 + (double)steps * akima->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( akima->cf(x).imag(), (f0 + (double)steps * akima->cdf(x0)).imag() );
	steps++;
	x += stepsize;
	EXPECT_DOUBLE_EQ( steffen->cf(x).real(), (f0 + (double)steps * steffen->cdf(x0)).real() );
	EXPECT_DOUBLE_EQ( steffen->cf(x).imag(), (f0 + (double)steps * steffen->cdf(x0)).imag() );


	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectComplexFirstDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	double x0 = linear->get_high_interp_limit();
	double step = 2.0;
	double xtest = x0 + step;
	EXPECT_DOUBLE_EQ( linear->df(xtest), linear->df(x0) );

	x0 = polynomial->get_high_interp_limit();
	xtest = x0 + step;
	EXPECT_DOUBLE_EQ( polynomial->cdf(xtest).real(), polynomial->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( cubic->cdf(xtest+1.0).real(), cubic->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( akima->cdf(xtest+1.5).real(), akima->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( steffen->cdf(xtest+2.5).real(), steffen->cdf(x0).real() );
	EXPECT_DOUBLE_EQ( polynomial->cdf(xtest).imag(), polynomial->cdf(x0).imag() );
	EXPECT_DOUBLE_EQ( cubic->cdf(xtest+1.0).imag(), cubic->cdf(x0).imag() );
	EXPECT_DOUBLE_EQ( akima->cdf(xtest+1.5).imag(), akima->cdf(x0).imag() );
	EXPECT_DOUBLE_EQ( steffen->cdf(xtest+2.5).imag(), steffen->cdf(x0).imag() );

	// test periodic in a minute
}

TEST_F(GSLInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectComplexSecondDerivative) {
	// Linear:
	// f(1) = 2.0
	// df/dx|x=1 = 2.0
	double x0 = linear->get_high_interp_limit();
	double step = 2.0;
	double xtest = x0 + step;
	EXPECT_DOUBLE_EQ( linear->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( linear->cdf(2,xtest).imag(), 0.0 );

	x0 = polynomial->get_high_interp_limit();
	xtest = x0 + step;
	EXPECT_DOUBLE_EQ( polynomial->cdf(2,xtest).real(), 0.0 );
	EXPECT_DOUBLE_EQ( cubic->cdf(2,xtest+1.0).real(), 0.0 );
	EXPECT_DOUBLE_EQ( akima->cdf(2,xtest+3.0).real(), 0.0 );
	EXPECT_DOUBLE_EQ( steffen->cdf(2,xtest+5.24).real(), 0.0 );
	EXPECT_DOUBLE_EQ( polynomial->cdf(2,xtest).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( cubic->cdf(2,xtest+1.0).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( akima->cdf(2,xtest+3.0).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( steffen->cdf(2,xtest+5.24).imag(), 0.0 );

	// test periodic in a minute
}






TEST_F(GSLInterpolator1DTest,LinearIdentifierIsCorrect) {
	EXPECT_EQ( linear->identifier(), Interpolator1D::as_string( linear->type() ) );

}

TEST_F(GSLInterpolator1DTest,PolynomialIdentifierIsCorrect) {
	EXPECT_EQ( polynomial->identifier(), Interpolator1D::as_string( polynomial->type() ) );
}

TEST_F(GSLInterpolator1DTest,CubicIdentifierIsCorrect) {
	EXPECT_EQ( cubic->identifier(), Interpolator1D::as_string( cubic->type() ) );

}

TEST_F(GSLInterpolator1DTest,AkimaIdentifierIsCorrect) {
	EXPECT_EQ( akima->identifier(), Interpolator1D::as_string( akima->type() ) );

}

TEST_F(GSLInterpolator1DTest,SteffenIdentifierIsCorrect) {
	EXPECT_EQ( steffen->identifier(), Interpolator1D::as_string( steffen->type() ) );
}


TEST_F(GSLInterpolator1DTest,InterpolatedValuesAreCorrectForLinear) {
	for (size_t i = 0; i < 4; i++) {
		EXPECT_NEAR( linear->f(l_interp_x[i]), l_interp_r[i], l_tolerance );
		EXPECT_NEAR( linear->cf(l_interp_x[i]).real(), l_interp_r[i], l_tolerance );
		EXPECT_NEAR( linear->cf(l_interp_x[i]).imag(), l_interp_i[i], l_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,InterpolatedValuesAreCorrectForPolynomial) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_NEAR( polynomial->f(n_interp_xvals[i]), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( polynomial->cf(n_interp_xvals[i]).real(), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( polynomial->cf(n_interp_xvals[i]).imag(), n_interp_ivals[i],
				fabs(n_interp_ivals[i])*n_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,InterpolatedValuesAreCorrectForCubic) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_NEAR( cubic->f(n_interp_xvals[i]), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( cubic->cf(n_interp_xvals[i]).real(), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( cubic->cf(n_interp_xvals[i]).imag(), n_interp_ivals[i],
				fabs(n_interp_ivals[i])*n_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,InterpolatedValuesAreCorrectForAkima) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_NEAR( akima->f(n_interp_xvals[i]), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( akima->cf(n_interp_xvals[i]).real(), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( akima->cf(n_interp_xvals[i]).imag(), n_interp_ivals[i],
				fabs(n_interp_ivals[i])*n_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,InterpolatedValuesAreCorrectForSteffen) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_NEAR( steffen->f(n_interp_xvals[i]), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( steffen->cf(n_interp_xvals[i]).real(), n_interp_rvals[i],
				fabs(n_interp_rvals[i])*n_tolerance );
		EXPECT_NEAR( steffen->cf(n_interp_xvals[i]).imag(), n_interp_ivals[i],
				fabs(n_interp_ivals[i])*n_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,FirstDerivativesAreCorrectForLinear) {
	for (size_t i = 0; i < 4; i++) {
		EXPECT_NEAR( linear->df(l_derivs_x[i]), l_derivs_r[i], l_tolerance );
		EXPECT_NEAR( linear->cdf(l_derivs_x[i]).real(), l_derivs_r[i], l_tolerance );
		EXPECT_NEAR( linear->cdf(l_derivs_x[i]).imag(), l_derivs_i[i], l_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,FirstDerivativesAreCorrectForPolynomial) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_NEAR( polynomial->df(n_deriv1_xvals[i]), n_deriv1_rvals[i],
				fabs(n_deriv1_rvals[i])*n_tolerance );
		EXPECT_NEAR( polynomial->cdf(n_deriv1_xvals[i]).real(), n_deriv1_rvals[i],
				fabs(n_deriv1_rvals[i])*n_tolerance );
		EXPECT_NEAR( polynomial->cdf(n_deriv1_xvals[i]).imag(), n_deriv1_ivals[i],
				fabs(n_deriv1_ivals[i])*n_tolerance );
	}
}

TEST_F(GSLInterpolator1DTest,FirstDerivativesHaveCorrectSignForCubic) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_GT( cubic->df(n_deriv1_xvals[i]), 0.0 );
		EXPECT_GT( cubic->cdf(n_deriv1_xvals[i]).real(), 0.0 );
		EXPECT_GT( cubic->cdf(n_deriv1_xvals[i]).imag(), 0.0 );
	}
}

TEST_F(GSLInterpolator1DTest,FirstDerivativesHaveCorrectSignForAkima) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_GT( akima->df(n_deriv1_xvals[i]), 0.0 );
		EXPECT_GT( akima->cdf(n_deriv1_xvals[i]).real(), 0.0 );
		EXPECT_GT( akima->cdf(n_deriv1_xvals[i]).imag(), 0.0 );
	}
}

TEST_F(GSLInterpolator1DTest,FirstDerivativesHaveCorrectSignForSteffen) {
	for (size_t i = 0; i < 8; i++) {
		EXPECT_GT( steffen->df(n_deriv1_xvals[i]), 0.0 );
		EXPECT_GT( steffen->cdf(n_deriv1_xvals[i]).real(), 0.0 );
		EXPECT_GT( steffen->cdf(n_deriv1_xvals[i]).imag(), 0.0 );
	}
}


TEST_F(GSLInterpolator1DTest,IsExtrapolatingReturnsTrueByDefault) {
	EXPECT_TRUE( linear->is_extrapolating() );
	EXPECT_TRUE( polynomial->is_extrapolating() );
	EXPECT_TRUE( cubic->is_extrapolating() );
	EXPECT_TRUE( akima->is_extrapolating() );
	EXPECT_TRUE( steffen->is_extrapolating() );
}

TEST_F(GSLInterpolator1DTest,IsExtrapolatingReturnsFalseIfSet) {
	linear->is_extrapolating( false );
	polynomial->is_extrapolating( false );
	cubic->is_extrapolating( false );
	akima->is_extrapolating( false );
	steffen->is_extrapolating( false );
	EXPECT_FALSE( linear->is_extrapolating() );
	EXPECT_FALSE( polynomial->is_extrapolating() );
	EXPECT_FALSE( cubic->is_extrapolating() );
	EXPECT_FALSE( akima->is_extrapolating() );
	EXPECT_FALSE( steffen->is_extrapolating() );
}

TEST_F(GSLInterpolator1DTest,ExtrapolationThrowsOutOfRangeIfExtrapolationDisabled) {
	Interpolator1D *interp1 = linear;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );

	interp1 = polynomial;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
	interp1 = cubic;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
	interp1 = akima;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
	interp1 = steffen;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
}


TEST_F(GSLInterpolator1DTest,ReturnsCorrectMaxDerivativeOrder) {
	EXPECT_EQ(cubic->max_derivative(),2);
}

// Tests that don't use the fixture
//TEST(SuiteName,TestName2) {
//
//}

/*
Example tests:
Integer Equality		EXPECT_EQ( v2.size(), 10 );
Floating Point Equality		EXPECT_DOUBLE_EQ( 40.0, 40.0 );
Vector contents			EXPECT_THAT( v2, ElementsAre(
					DoubleEq(1.0),
					DoubleEq(2.0)) );
Exception Thrown Properly	EXPECT_THROW( {result = command();},
					std::out_of_range );
*/
