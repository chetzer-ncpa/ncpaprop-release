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
class GSLInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		linear = Interpolator1D::build(interpolator1d_t::GSL_1D_LINEAR);
		polynomial = Interpolator1D::build(interpolator1d_t::GSL_1D_POLYNOMIAL);
		cubic = Interpolator1D::build(interpolator1d_t::GSL_1D_CUBIC);
		cubic_periodic = Interpolator1D::build(interpolator1d_t::GSL_1D_CUBIC_PERIODIC);
		akima = Interpolator1D::build(interpolator1d_t::GSL_1D_AKIMA);
		akima_periodic = Interpolator1D::build(interpolator1d_t::GSL_1D_AKIMA_PERIODIC);
		steffen = Interpolator1D::build(interpolator1d_t::GSL_1D_STEFFEN);

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
				n_deriv1_xvals[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_deriv1_rvals[i-1] = evalpoly( 3, coeffs1d, n_deriv1_xvals[i-1] );
				n_deriv1_ivals[i-1] = 0.5 * n_deriv1_rvals[i-1];
				n_deriv2_xvals[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_deriv2_rvals[i-1] = evalpoly( 2, coeffs2d, n_deriv2_xvals[i-1] );
				n_deriv2_ivals[i-1] = 0.5 * n_deriv2_rvals[i-1];
			}

			// test function for nonlinear periodic interpolators is
			// y = cos((x-1)*PI/4.0) = cos( x*PI/4 - PI/4 )
			// dy/dx = -(PI/4)sin( x*PI/4 - PI/4 )
			// d2y/dx2 = -(PI**2/16)cos( x*PI/4 - PI/4 )
			n_rvals_p[i] = std::cos( (n_xvals[i] - 1)*phase );
			n_ivals_p[i] = 0.5 * n_rvals_p[i];
			if (i > 0) {
				n_deriv1_xvals_p[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_deriv1_rvals_p[i-1] = -phase * std::sin( (n_deriv1_xvals_p[i-1] - 1)*phase );
				n_deriv1_ivals_p[i-1] = 0.5 * n_deriv1_rvals_p[i-1];
				n_deriv2_xvals_p[i-1] = 0.5 * (n_xvals[i] + n_xvals[i-1]);
				n_deriv2_rvals_p[i-1] = -phase * phase * std::cos( (n_deriv1_xvals_p[i-1] - 1)*phase );
				n_deriv2_ivals_p[i-1] = 0.5 * n_deriv2_rvals_p[i-1];
			}
		}
		linear->set( 5, l_xvals, l_rvals, l_ivals );
		polynomial->set( 9, n_xvals, n_rvals, n_ivals );
		cubic->set( 9, n_xvals, n_rvals, n_ivals );
		akima->set( 9, n_xvals, n_rvals, n_ivals );
		steffen->set( 9, n_xvals, n_rvals, n_ivals );
		cubic_periodic->set( 9, n_xvals, n_rvals_p, n_ivals_p );
		akima_periodic->set( 9, n_xvals, n_rvals_p, n_ivals_p );

		double2complex( 5, l_rvals, l_ivals, l_cvals );
	}

	// If there's any cleanup other than normal destructors
//	void TearDown() override {
//		delete interp1;
//		delete c_interp1;
//	}

	// class members here
	Interpolator1D *linear, *polynomial, *cubic, *cubic_periodic,
					*akima, *akima_periodic, *steffen;
	GSLInterpolator1D *gPtr1, *gPtr2;

	// vectors for testing linear interpolators
	double l_xvals[5] = { 1, 2, 3, 4, 5 };
	double l_rvals[5] = { 2, 4, 8, 16, 32 };
	double l_ivals[5] = { 2, 1, 0, -2, -4 };
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
	double n_rvals[9], n_rvals_p[9], n_ivals[9], n_ivals_p[9];
	complex<double> n_cvals[9], n_cvals_p[9];
	double n_deriv1_xvals[8], n_deriv1_rvals[8], n_deriv1_ivals[8];
	double n_deriv2_xvals[8], n_deriv2_rvals[8], n_deriv2_ivals[8];
	double n_deriv1_xvals_p[8], n_deriv1_rvals_p[8], n_deriv1_ivals_p[8];
	double n_deriv2_xvals_p[8], n_deriv2_rvals_p[8], n_deriv2_ivals_p[8];


	double zvals5[5] = { 0, 0, 0, 0, 0 };
	double zvals9[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
};


// Tests that use the fixture
TEST_F(GSLInterpolator1DTest,CopyConstructorCreatesCopy) {
	complex<double> testval = linear->eval_cf( 2.1 );
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	GSLInterpolator1D linear2( *gPtr1 );

	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_r_->x, gPtr1->spline_r_->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_r_->y, gPtr1->spline_r_->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_i_->x, gPtr1->spline_i_->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_i_->y, gPtr1->spline_i_->y );
	EXPECT_EQ( linear2.accel_r_->cache, gPtr1->accel_r_->cache );
	EXPECT_EQ( linear2.accel_r_->miss_count, gPtr1->accel_r_->miss_count );
	EXPECT_EQ( linear2.accel_r_->hit_count, gPtr1->accel_r_->hit_count );
	EXPECT_EQ( linear2.accel_i_->cache, gPtr1->accel_i_->cache );
	EXPECT_EQ( linear2.accel_i_->miss_count, gPtr1->accel_i_->miss_count );
	EXPECT_EQ( linear2.accel_i_->hit_count, gPtr1->accel_i_->hit_count );
	EXPECT_TRUE( linear2.is_ready() );
}

TEST_F(GSLInterpolator1DTest,SwapWorksProperly) {
	double testval = steffen->eval_f( 2.5 );
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	gPtr2 = static_cast<GSLInterpolator1D *>( steffen );
	GSLInterpolator1D linear2( *gPtr1 ),
					  steffen2( *gPtr2 );
	EXPECT_EQ( linear2.accel_r_->cache, gPtr1->accel_r_->cache );
	EXPECT_EQ( steffen2.accel_r_->cache, gPtr2->accel_r_->cache );
	EXPECT_EQ( linear2.accel_r_->miss_count, gPtr1->accel_r_->miss_count );
	EXPECT_EQ( steffen2.accel_r_->miss_count, gPtr2->accel_r_->miss_count );
	EXPECT_EQ( linear2.accel_r_->hit_count, gPtr1->accel_r_->hit_count );
	EXPECT_EQ( steffen2.accel_r_->hit_count, gPtr2->accel_r_->hit_count );

	swap( linear2, steffen2 );

	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_r_->x, gPtr2->spline_r_->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_r_->y, gPtr2->spline_r_->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_i_->x, gPtr2->spline_i_->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, linear2.spline_i_->y, gPtr2->spline_i_->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.spline_r_->x, gPtr1->spline_r_->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.spline_r_->y, gPtr1->spline_r_->y );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.spline_i_->x, gPtr1->spline_i_->x );
	EXPECT_DOUBLE_ARRAY_EQ( 5, steffen2.spline_i_->y, gPtr1->spline_i_->y );
	EXPECT_EQ( linear2.identifier(), gPtr2->identifier() );
	EXPECT_EQ( steffen2.identifier(), gPtr1->identifier() );
	EXPECT_EQ( linear2.accel_r_->cache, gPtr2->accel_r_->cache );
	EXPECT_EQ( steffen2.accel_r_->cache, gPtr1->accel_r_->cache );
	EXPECT_EQ( linear2.accel_r_->miss_count, gPtr2->accel_r_->miss_count );
	EXPECT_EQ( steffen2.accel_r_->miss_count, gPtr1->accel_r_->miss_count );
	EXPECT_EQ( linear2.accel_r_->hit_count, gPtr2->accel_r_->hit_count );
	EXPECT_EQ( steffen2.accel_r_->hit_count, gPtr1->accel_r_->hit_count );
}

TEST_F(GSLInterpolator1DTest,RealSetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_r_->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_r_->y, l_rvals );
	gPtr1->set( 9, n_xvals, n_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_r_->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_r_->y, n_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_i_->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_i_->y, zvals9 );

}

TEST_F(GSLInterpolator1DTest,RealAndImaginarySetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLInterpolator1D *>( linear );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_r_->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_r_->y, l_rvals );
	gPtr1->set( 9, n_xvals, n_rvals, n_ivals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_r_->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_r_->y, n_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_i_->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_i_->y, n_ivals );
}

TEST_F(GSLInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	gPtr1 = static_cast<GSLInterpolator1D *>( steffen );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_r_->x, n_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 9, gPtr1->spline_r_->y, n_rvals );
	gPtr1->set( 5, l_xvals, l_cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_r_->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_r_->y, l_rvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_i_->x, l_xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, gPtr1->spline_i_->y, l_ivals );
}

TEST_F(GSLInterpolator1DTest,InitSetsPointersToNull) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	gPtr1->init();
	EXPECT_EQ( gPtr1->spline_r_, nullptr );
	EXPECT_EQ( gPtr1->spline_i_, nullptr );
	EXPECT_EQ( gPtr1->accel_r_, nullptr );
	EXPECT_EQ( gPtr1->accel_i_, nullptr );
}

TEST_F(GSLInterpolator1DTest,AllocateSetsPointersToValid) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	gPtr1->init();
	gPtr1->allocate(15);
	EXPECT_NE( gPtr1->spline_r_, nullptr );
	EXPECT_NE( gPtr1->spline_i_, nullptr );
	EXPECT_NE( gPtr1->accel_r_, nullptr );
	EXPECT_NE( gPtr1->accel_i_, nullptr );
}

TEST_F(GSLInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
	gPtr1 = static_cast<GSLInterpolator1D *>( cubic );
	gPtr1->init();
}

TEST_F(GSLInterpolator1DTest,FreeSetsPointersToNull) {

}

TEST_F(GSLInterpolator1DTest,IsReadyReturnsProperValue) {

}

TEST_F(GSLInterpolator1DTest,IdentifierIsCorrect) {

}

TEST_F(GSLInterpolator1DTest,RealInterpolatedValuesAreCorrect) {
}

TEST_F(GSLInterpolator1DTest,ComplexInterpolatedValuesAreCorrect) {

}

TEST_F(GSLInterpolator1DTest,RealInterpolatedFirstDerivativesAreCorrect) {

}

TEST_F(GSLInterpolator1DTest,ComplexInterpolatedFirstDerivativesAreCorrect) {

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
