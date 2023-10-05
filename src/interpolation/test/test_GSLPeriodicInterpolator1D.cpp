#include "NCPAInterpolation.h"
#include "NCPACommon.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <complex>

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
class GSLPeriodicInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		double phase = M_PI / 40.0;
		for (size_t i = 0; i < 81; i++) {
			n_xvals_p[i] = 1.0 + ((double)i) / 10.0;

			// test function for nonlinear periodic interpolators is
			// y = cos((x-1)*PI/40.0) = cos( x*PI/40 - PI/40 )
			// dy/dx = -(PI/40)sin( x*PI/40 - PI/40 )
			// d2y/dx2 = -(PI**2/1600)cos( x*PI/40 - PI/40 )
			n_rvals_p[i] = std::cos( (n_rvals_p[i] - 1)*phase );
			n_ivals_p[i] = 0.5 * n_rvals_p[i];
			n_cvals_p[i].real(n_rvals_p[i]);
			n_cvals_p[i].imag(n_ivals_p[i]);
			if (i > 0) {
				n_interp_xvals_p[i-1] = 0.5 * (n_rvals_p[i] + n_rvals_p[i-1]);
				n_interp_rvals_p[i-1] = std::cos( (n_interp_xvals_p[i-1] - 1)*phase );
				n_interp_ivals_p[i-1] = 0.5 * n_interp_rvals_p[i-1];
				n_deriv1_xvals_p[i-1] = 0.5 * (n_xvals_p[i] + n_xvals_p[i-1]);
				n_deriv1_rvals_p[i-1] = -phase * std::sin( (n_deriv1_xvals_p[i-1] - 1)*phase );
				n_deriv1_ivals_p[i-1] = 0.5 * n_deriv1_rvals_p[i-1];
				n_deriv2_xvals_p[i-1] = 0.5 * (n_xvals_p[i] + n_xvals_p[i-1]);
				n_deriv2_rvals_p[i-1] = -phase * phase * std::cos( (n_deriv1_xvals_p[i-1] - 1)*phase );
				n_deriv2_ivals_p[i-1] = 0.5 * n_deriv2_rvals_p[i-1];
			}
		}

		for (size_t i = 0; i < 9; i++) {
			n_xvals_p50[i] = n_xvals_p[i+50];
			n_rvals_p50[i] = n_rvals_p[i+50];
			n_ivals_p50[i] = n_ivals_p[i+50];
		}

		cubic_periodic = Interpolator1D::build(interpolator1d_t::GSL_1D_CUBIC_PERIODIC);
		akima_periodic = Interpolator1D::build(interpolator1d_t::GSL_1D_AKIMA_PERIODIC);

		cubic_periodic->set( 81, n_xvals_p, n_rvals_p, n_ivals_p )->ready();
		akima_periodic->set( 81, n_xvals_p, n_rvals_p, n_ivals_p )->ready();

	}

	// If there's any cleanup other than normal destructors
//	void TearDown() override {
//		delete interp1;
//		delete c_interp1;
//	}

	// class members here
	Interpolator1D *cubic_periodic, *akima_periodic;
	GSLPeriodicInterpolator1D *gPtr1, *gPtr2;
	double l_tolerance = 1e-3;
	double n_tolerance = 0.15;

	// vectors for testing nonlinear interpolators
	double n_xvals_p[81];
//	double coeffs[4] = { 2.0, -1.5, 1.25, 2.0 };
//	double coeffs1d[3] = { -1.5, 2.5, 6.0 };
//	double coeffs2d[2] = { 2.5, 12.0 };

	// to calculate
	double n_rvals[81], n_rvals_p[81], n_ivals[81], n_ivals_p[81];
	complex<double> n_cvals[81], n_cvals_p[81];
	double n_interp_xvals_p[80], n_interp_rvals_p[80], n_interp_ivals_p[80];
	double n_deriv1_xvals_p[80], n_deriv1_rvals_p[80], n_deriv1_ivals_p[80];
	double n_deriv2_xvals_p[80], n_deriv2_rvals_p[80], n_deriv2_ivals_p[80];
	double n_xvals_p50[9], n_rvals_p50[9], n_ivals_p50[9];


	double zvals5[5] = { 0, 0, 0, 0, 0 };
	double zvals9[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double zvals81[81];
};


// Tests that use the fixture
TEST_F(GSLPeriodicInterpolator1DTest,CopyConstructorCreatesCopy) {
	complex<double> testval = cubic_periodic->cf( 2.1 );
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	GSLPeriodicInterpolator1D cubic_periodic2( *gPtr1 );

	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic_periodic2.get_real_spline()->x, gPtr1->get_real_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic_periodic2.get_real_spline()->y, gPtr1->get_real_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic_periodic2.get_imag_spline()->x, gPtr1->get_imag_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic_periodic2.get_imag_spline()->y, gPtr1->get_imag_spline()->y );
	EXPECT_EQ( cubic_periodic2.get_real_accel()->cache, gPtr1->get_real_accel()->cache );
	EXPECT_EQ( cubic_periodic2.get_real_accel()->miss_count, gPtr1->get_real_accel()->miss_count );
	EXPECT_EQ( cubic_periodic2.get_real_accel()->hit_count, gPtr1->get_real_accel()->hit_count );
	EXPECT_EQ( cubic_periodic2.get_imag_accel()->cache, gPtr1->get_imag_accel()->cache );
	EXPECT_EQ( cubic_periodic2.get_imag_accel()->miss_count, gPtr1->get_imag_accel()->miss_count );
	EXPECT_EQ( cubic_periodic2.get_imag_accel()->hit_count, gPtr1->get_imag_accel()->hit_count );
	EXPECT_TRUE( cubic_periodic2.is_ready() );
}

TEST_F(GSLPeriodicInterpolator1DTest,SwapWorksProperly) {
	double testval = cubic_periodic->f( 2.5 );
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	gPtr2 = static_cast<GSLPeriodicInterpolator1D *>( akima_periodic );
	GSLPeriodicInterpolator1D cubic2( *gPtr1 ),
					  akima2( *gPtr2 );
	EXPECT_EQ( cubic2.get_real_accel()->cache, gPtr1->get_real_accel()->cache );
	EXPECT_EQ( akima2.get_real_accel()->cache, gPtr2->get_real_accel()->cache );
	EXPECT_EQ( cubic2.get_real_accel()->miss_count, gPtr1->get_real_accel()->miss_count );
	EXPECT_EQ( akima2.get_real_accel()->miss_count, gPtr2->get_real_accel()->miss_count );
	EXPECT_EQ( cubic2.get_real_accel()->hit_count, gPtr1->get_real_accel()->hit_count );
	EXPECT_EQ( akima2.get_real_accel()->hit_count, gPtr2->get_real_accel()->hit_count );

	swap( cubic2, akima2 );

	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic2.get_real_spline()->x, gPtr2->get_real_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic2.get_real_spline()->y, gPtr2->get_real_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic2.get_imag_spline()->x, gPtr2->get_imag_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 81, cubic2.get_imag_spline()->y, gPtr2->get_imag_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 81, akima2.get_real_spline()->x, gPtr1->get_real_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 81, akima2.get_real_spline()->y, gPtr1->get_real_spline()->y );
	EXPECT_DOUBLE_ARRAY_EQ( 81, akima2.get_imag_spline()->x, gPtr1->get_imag_spline()->x );
	EXPECT_DOUBLE_ARRAY_EQ( 81, akima2.get_imag_spline()->y, gPtr1->get_imag_spline()->y );
	EXPECT_EQ( cubic2.identifier(), gPtr2->identifier() );
	EXPECT_EQ( akima2.identifier(), gPtr1->identifier() );
	EXPECT_EQ( cubic2.get_real_accel()->cache, gPtr2->get_real_accel()->cache );
	EXPECT_EQ( akima2.get_real_accel()->cache, gPtr1->get_real_accel()->cache );
	EXPECT_EQ( cubic2.get_real_accel()->miss_count, gPtr2->get_real_accel()->miss_count );
	EXPECT_EQ( akima2.get_real_accel()->miss_count, gPtr1->get_real_accel()->miss_count );
	EXPECT_EQ( cubic2.get_real_accel()->hit_count, gPtr2->get_real_accel()->hit_count );
	EXPECT_EQ( akima2.get_real_accel()->hit_count, gPtr1->get_real_accel()->hit_count );
}

TEST_F(GSLPeriodicInterpolator1DTest,CloneWorksProperly) {
	Interpolator1D *clone = cubic_periodic->clone();
	cubic_periodic->free();
	for (size_t i = 0; i < 81; i++) {
		EXPECT_DOUBLE_EQ( clone->f(n_xvals_p[i]), n_rvals_p[i] );
	}
	delete clone;
	clone = akima_periodic->clone();
	akima_periodic->free();
	for (size_t i = 0; i < 81; i++) {
		EXPECT_DOUBLE_EQ( clone->f(n_xvals_p[i]), n_rvals_p[i] );
	}
	delete clone;
}

TEST_F(GSLPeriodicInterpolator1DTest,RealSetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	EXPECT_DOUBLE_ARRAY_EQ( 81, gPtr1->get_real_spline()->x, n_xvals_p );
	EXPECT_DOUBLE_ARRAY_EQ( 81, gPtr1->get_real_spline()->y, n_rvals_p );
	gPtr1->set( 9, n_xvals_p50, n_rvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->x, n_xvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->y, n_rvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->x, n_xvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->y, zvals9 );

}

TEST_F(GSLPeriodicInterpolator1DTest,RealAndImaginarySetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	EXPECT_DOUBLE_ARRAY_EQ( 81, gPtr1->get_real_spline()->x, n_xvals_p );
	EXPECT_DOUBLE_ARRAY_EQ( 81, gPtr1->get_real_spline()->y, n_rvals_p );
	gPtr1->set( 9, n_xvals_p50, n_rvals_p50, n_ivals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->x, n_xvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_real_spline()->y, n_rvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->x, n_xvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->get_imag_spline()->y, n_ivals_p50 );
}

TEST_F(GSLPeriodicInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	EXPECT_DOUBLE_ARRAY_EQ( 81, gPtr1->get_real_spline()->x, n_xvals_p );
	EXPECT_DOUBLE_ARRAY_EQ( 81, gPtr1->get_real_spline()->y, n_rvals_p );
	gPtr1->set( 5, n_xvals_p50, n_cvals_p+50 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->x, n_xvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_real_spline()->y, n_rvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_imag_spline()->x, n_xvals_p50 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->get_imag_spline()->y, n_ivals_p50 );
}

TEST_F(GSLPeriodicInterpolator1DTest,InitSetsPointersToNull) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	gPtr1->init();
	EXPECT_EQ( gPtr1->get_real_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_real_accel(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_accel(), nullptr );
}

TEST_F(GSLPeriodicInterpolator1DTest,AllocateSetsPointersToValid) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	gPtr1->init();
	gPtr1->allocate(15);
	EXPECT_NE( gPtr1->get_real_spline(), nullptr );
	EXPECT_NE( gPtr1->get_imag_spline(), nullptr );
	EXPECT_NE( gPtr1->get_real_accel(), nullptr );
	EXPECT_NE( gPtr1->get_imag_accel(), nullptr );
}

TEST_F(GSLPeriodicInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	gPtr1->init();
}

TEST_F(GSLPeriodicInterpolator1DTest,FreeSetsPointersToNull) {
	gPtr1 = static_cast<GSLPeriodicInterpolator1D *>( cubic_periodic );
	cubic_periodic->free();
	EXPECT_EQ( gPtr1->get_real_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_spline(), nullptr );
	EXPECT_EQ( gPtr1->get_real_accel(), nullptr );
	EXPECT_EQ( gPtr1->get_imag_accel(), nullptr );
}

TEST_F(GSLPeriodicInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( cubic_periodic->is_ready() );
}

TEST_F(GSLPeriodicInterpolator1DTest,GetInterpLimitsReturnsCorrectLimits) {
	double x1, x2;
	cubic_periodic->get_interp_limits(x1, x2);
	EXPECT_DOUBLE_EQ( x1, 1.0 );
	EXPECT_DOUBLE_EQ( x2, 9.0 );
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectRealValue) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->f( testx - xrange ), cubic_periodic->f( testx ));
	EXPECT_DOUBLE_EQ( akima_periodic->f( testx - xrange ), akima_periodic->f( testx ));
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectRealFirstDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->df( testx - xrange ), cubic_periodic->df( testx ));
	EXPECT_DOUBLE_EQ( akima_periodic->df( testx - xrange ), akima_periodic->df( testx ));
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectRealSecondDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->df( 2, testx - xrange ), cubic_periodic->df( 2, testx ));
	EXPECT_DOUBLE_EQ( akima_periodic->df( 2, testx - xrange ), akima_periodic->df( 2, testx ));
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectRealValue) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->f( testx + xrange ), cubic_periodic->f( testx ));
	EXPECT_DOUBLE_EQ( akima_periodic->f( testx + xrange ), akima_periodic->f( testx ));
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectRealFirstDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->df( testx + xrange ), cubic_periodic->df( testx ));
	EXPECT_DOUBLE_EQ( akima_periodic->df( testx + xrange ), akima_periodic->df( testx ));
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectRealSecondDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->df( 2, testx + xrange ), cubic_periodic->df( 2, testx ));
	EXPECT_DOUBLE_EQ( akima_periodic->df( 2, testx + xrange ), akima_periodic->df( 2, testx ));
}


TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectComplexValue) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->cf( testx - xrange ).real(),
			cubic_periodic->cf( testx ).real());
	EXPECT_DOUBLE_EQ( cubic_periodic->cf( testx - xrange ).imag(),
			cubic_periodic->cf( testx ).imag());
	EXPECT_DOUBLE_EQ( akima_periodic->cf( testx - xrange ).real(),
			akima_periodic->cf( testx ).real());
	EXPECT_DOUBLE_EQ( akima_periodic->cf( testx - xrange ).imag(),
			akima_periodic->cf( testx ).imag());
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectComplexFirstDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( testx - xrange ).real(),
			cubic_periodic->cdf( testx ).real());
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( testx - xrange ).imag(),
			cubic_periodic->cdf( testx ).imag());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( testx - xrange ).real(),
			akima_periodic->cdf( testx ).real());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( testx - xrange ).imag(),
			akima_periodic->cdf( testx ).imag());
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnLowEndReturnsCorrectComplexSecondDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( 2, testx - xrange ).real(),
			cubic_periodic->cdf( 2, testx ).real());
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( 2, testx - xrange ).imag(),
			cubic_periodic->cdf( 2, testx ).imag());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( 2, testx - xrange ).real(),
			akima_periodic->cdf( 2, testx ).real());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( 2, testx - xrange ).imag(),
			akima_periodic->cdf( 2, testx ).imag());
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectComplexValue) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->cf( testx + xrange ).real(),
			cubic_periodic->cf( testx ).real());
	EXPECT_DOUBLE_EQ( cubic_periodic->cf( testx + xrange ).imag(),
			cubic_periodic->cf( testx ).imag());
	EXPECT_DOUBLE_EQ( akima_periodic->cf( testx + xrange ).real(),
			akima_periodic->cf( testx ).real());
	EXPECT_DOUBLE_EQ( akima_periodic->cf( testx + xrange ).imag(),
			akima_periodic->cf( testx ).imag());
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectComplexFirstDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( testx + xrange ).real(),
			cubic_periodic->cdf( testx ).real());
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( testx + xrange ).imag(),
			cubic_periodic->cdf( testx ).imag());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( testx + xrange ).real(),
			akima_periodic->cdf( testx ).real());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( testx + xrange ).imag(),
			akima_periodic->cdf( testx ).imag());
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationOnHighEndReturnsCorrectComplexSecondDerivative) {
	double xrange = cubic_periodic->get_high_interp_limit()
			- cubic_periodic->get_low_interp_limit();
	double testx = 3.0;
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( 2, testx + xrange ).real(),
			cubic_periodic->cdf( 2, testx ).real());
	EXPECT_DOUBLE_EQ( cubic_periodic->cdf( 2, testx + xrange ).imag(),
			cubic_periodic->cdf( 2, testx ).imag());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( 2, testx + xrange ).real(),
			akima_periodic->cdf( 2, testx ).real());
	EXPECT_DOUBLE_EQ( akima_periodic->cdf( 2, testx + xrange ).imag(),
			akima_periodic->cdf( 2, testx ).imag());
}



TEST_F(GSLPeriodicInterpolator1DTest,CubicPeriodicIdentifierIsCorrect) {
	EXPECT_EQ( cubic_periodic->identifier(), Interpolator1D::as_string( cubic_periodic->type() ) );

}


TEST_F(GSLPeriodicInterpolator1DTest,AkimaPeriodicIdentifierIsCorrect) {
	EXPECT_EQ( akima_periodic->identifier(), Interpolator1D::as_string( akima_periodic->type() ) );

}

TEST_F(GSLPeriodicInterpolator1DTest,InterpolatedValuesAreCorrectForCubicPeriodic) {
	for (size_t i = 0; i < 80; i++) {
		EXPECT_NEAR( cubic_periodic->f(n_interp_xvals_p[i]), n_interp_rvals_p[i],
				fabs(n_interp_rvals_p[i])*n_tolerance );
		EXPECT_NEAR( cubic_periodic->cf(n_interp_xvals_p[i]).real(), n_interp_rvals_p[i],
				fabs(n_interp_rvals_p[i])*n_tolerance );
		EXPECT_NEAR( cubic_periodic->cf(n_interp_xvals_p[i]).imag(), n_interp_ivals_p[i],
				fabs(n_interp_ivals_p[i])*n_tolerance );
	}
}

TEST_F(GSLPeriodicInterpolator1DTest,InterpolatedValuesAreCorrectForAkimaPeriodic) {
	for (size_t i = 0; i < 80; i++) {
		EXPECT_NEAR( akima_periodic->f(n_interp_xvals_p[i]), n_interp_rvals_p[i],
				fabs(n_interp_rvals_p[i])*n_tolerance );
		EXPECT_NEAR( akima_periodic->cf(n_interp_xvals_p[i]).real(), n_interp_rvals_p[i],
				fabs(n_interp_rvals_p[i])*n_tolerance );
		EXPECT_NEAR( akima_periodic->cf(n_interp_xvals_p[i]).imag(), n_interp_ivals_p[i],
				fabs(n_interp_ivals_p[i])*n_tolerance );
	}
}



TEST_F(GSLPeriodicInterpolator1DTest,FirstDerivativesCloseEnoughForCubicPeriodic) {
	for (size_t i = 0; i < 80; i++) {
		EXPECT_NEAR( cubic_periodic->df(n_deriv1_xvals_p[i]), n_deriv1_rvals_p[i],
				0.05 );
		EXPECT_NEAR( cubic_periodic->cdf(n_deriv1_xvals_p[i]).real(), n_deriv1_rvals_p[i],
				0.05 );
		EXPECT_NEAR( cubic_periodic->cdf(n_deriv1_xvals_p[i]).imag(), n_deriv1_ivals_p[i],
				0.05 );
	}
}

TEST_F(GSLPeriodicInterpolator1DTest,FirstDerivativesCloseEnoughForAkimaPeriodic) {
	for (size_t i = 0; i < 80; i++) {
		EXPECT_NEAR( akima_periodic->df(n_deriv1_xvals_p[i]), n_deriv1_rvals_p[i],
				0.05 );
		EXPECT_NEAR( akima_periodic->cdf(n_deriv1_xvals_p[i]).real(), n_deriv1_rvals_p[i],
				0.05 );
		EXPECT_NEAR( akima_periodic->cdf(n_deriv1_xvals_p[i]).imag(), n_deriv1_ivals_p[i],
				0.05 );
	}
}

TEST_F(GSLPeriodicInterpolator1DTest,IsExtrapolatingReturnsTrueByDefault) {
	EXPECT_TRUE( cubic_periodic->is_extrapolating() );
	EXPECT_TRUE( akima_periodic->is_extrapolating() );
}

TEST_F(GSLPeriodicInterpolator1DTest,IsExtrapolatingReturnsFalseIfSet) {
	cubic_periodic->is_extrapolating( false );
	akima_periodic->is_extrapolating( false );
	EXPECT_FALSE( cubic_periodic->is_extrapolating() );
	EXPECT_FALSE( akima_periodic->is_extrapolating() );
}

TEST_F(GSLPeriodicInterpolator1DTest,ExtrapolationThrowsOutOfRangeIfExtrapolationDisabled) {
	Interpolator1D *interp1 = cubic_periodic;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );

	interp1 = akima_periodic;
	interp1->is_extrapolating(false);
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );

}


TEST_F(GSLPeriodicInterpolator1DTest,ReturnsCorrectMaxDerivativeOrder) {
	EXPECT_EQ(cubic_periodic->max_derivative(),2);
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
