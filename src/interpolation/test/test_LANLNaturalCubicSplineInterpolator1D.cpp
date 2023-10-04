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
class LANLNaturalCubicSplineInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		// make cubic function to test with
		for (size_t i = 0; i < 5; i++) {
			yvals[i] = NCPA::evalpoly( 4, coeffs, xvals[i] );
			ivals[i] = 0.5 * yvals[i];
			cvals[i] = complex<double>(yvals[i],ivals[i]);

		}
		for (size_t i = 0; i < 4; i++) {
			double x = 0.5*(xvals[i] + xvals[i+1]);
			ymids[i] = NCPA::evalpoly( 4, coeffs, x );
			imids[i] = 0.5*ymids[i];
			cmids[i] = complex<double>(ymids[i],imids[i]);
			derivs1[i] = NCPA::evalpoly( 3, coeffs1d, x);
			iderivs1[i] = 0.5 * derivs1[i];
			derivs2[i] = NCPA::evalpoly( 2, coeffs2d, x );
			iderivs2[i] = 0.5 * derivs2[i];
			derivs3[i] = NCPA::evalpoly( 1, coeffs3d, x );
			iderivs3[i] = 0.5 * derivs3[i];

		}

		interp1 = Interpolator1D::build( interpolator1d_t::LANL_1D_NATURAL_CUBIC );
		interp1->set( 5, xvals, yvals )->ready();

		c_interp1 = Interpolator1D::build( interpolator1d_t::LANL_1D_NATURAL_CUBIC );
		c_interp1->set( 5, xvals, cvals )->ready();

		pInterp1 = static_cast<LANLNaturalCubicSplineInterpolator1D*>(interp1);
		pCInterp1 = static_cast<LANLNaturalCubicSplineInterpolator1D*>(c_interp1);
	}

	// If there's any cleanup other than normal destructors
	void TearDown() override {
		delete interp1;
		delete c_interp1;
	}

	// class members here
	Interpolator1D *interp1, *c_interp1;
	LANLNaturalCubicSplineInterpolator1D *pInterp1, *pCInterp1;
	double xvals[5] = { 1, 2, 3, 4, 5 };
	double yvals[5], ymids[4];
	double zvals[5] = { 0, 0, 0, 0, 0 };
	double ivals[5], imids[4];
	complex<double> cvals[5], cmids[4];
	size_t order = 3;
	double coeffs[4] = { 2.0, -1.5, 1.25, 2.0 };
	double coeffs1d[3] = { -1.5, 2.5, 6.0 };
	double coeffs2d[2] = { 2.5, 12.0 };
	double coeffs3d[1] = { 12.0 };
	double derivs1[4], derivs2[4], derivs3[4];
	double iderivs1[4], iderivs2[4], iderivs3[4];

	double tolerance = 0.1;
};


// Tests that use the fixture
TEST_F(LANLNaturalCubicSplineInterpolator1DTest,CopyConstructorCreatesCopy) {
	LANLNaturalCubicSplineInterpolator1D interp2( *pInterp1 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->x_vals,
			pInterp1->get_real_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->f_vals,
			pInterp1->get_real_spline()->f_vals );
	EXPECT_EQ( interp2.is_ready(), pInterp1->is_ready() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,SwapWorksProperly) {
	LANLNaturalCubicSplineInterpolator1D interp2( *pInterp1 ), c_interp2(*pCInterp1);
	swap( interp2, c_interp2 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_real_spline()->x_vals,
			pInterp1->get_real_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_real_spline()->f_vals,
			pInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->x_vals,
			pCInterp1->get_real_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->f_vals,
			pCInterp1->get_real_spline()->f_vals );

	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_imag_spline()->x_vals,
			pCInterp1->get_imag_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_imag_spline()->f_vals,
			pCInterp1->get_imag_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_imag_spline()->x_vals,
			pInterp1->get_imag_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_imag_spline()->f_vals,
			pInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,RealSetOverwritesPriorValues) {
	c_interp1->set( 5, xvals, xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, xvals, pCInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pCInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	interp1->set( 5, xvals, cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pCInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pCInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,SetXOverwritesValuesAndMakesUnready) {
	pInterp1->set_x( 5, yvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pInterp1->get_real_spline()->x_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,RealSetYOverwritesValuesAndMakesUnready) {
	pInterp1->set_y( 5, zvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pInterp1->get_real_spline()->f_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ComplexSetYOverwritesValuesAndMakesUnready) {
	pInterp1->set_y( 5, cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pInterp1->get_imag_spline()->f_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,InitClearsSplines) {
	interp1->init();
	EXPECT_EQ( pInterp1->get_real_spline()->length,0);
	EXPECT_EQ( pInterp1->get_imag_spline()->length,0);
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,InitSetsPointersToNull) {
	interp1->init();
	EXPECT_EQ( pInterp1->get_real_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_real_spline()->f_vals, nullptr);
	EXPECT_EQ( pInterp1->get_imag_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_imag_spline()->f_vals, nullptr);
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,AllocateSetsPointersToValid) {
	interp1->init();
	interp1->allocate( 15 );
	EXPECT_NE( pInterp1->get_real_spline()->x_vals, nullptr );
	EXPECT_NE( pInterp1->get_real_spline()->f_vals, nullptr);
	EXPECT_NE( pInterp1->get_imag_spline()->x_vals, nullptr );
	EXPECT_NE( pInterp1->get_imag_spline()->f_vals, nullptr);
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
	interp1->init();
	EXPECT_FALSE( interp1->is_ready() );
	interp1->allocate(5);
	EXPECT_FALSE( interp1->is_ready() );
	pInterp1->set_x( 5, yvals );
	EXPECT_FALSE( interp1->is_ready() );
	pInterp1->set_y( 5, zvals );
	EXPECT_FALSE( interp1->is_ready() );
	interp1->ready();
	EXPECT_TRUE( interp1->is_ready() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,FreeSetsPointersToNull) {
	interp1->free();
	EXPECT_EQ( pInterp1->get_real_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_real_spline()->f_vals, nullptr);
	EXPECT_EQ( pInterp1->get_imag_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_imag_spline()->f_vals, nullptr);
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( interp1->is_ready() );
	pInterp1->free();
	EXPECT_FALSE( interp1->is_ready() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,IdentifierIsCorrect) {
	EXPECT_EQ( interp1->identifier(), Interpolator1D::as_string( interp1->type() ));
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,RealInterpolatedValuesAreCorrect) {
	for (size_t i = 0; i < 4; i++) {
		double x = 0.5*(xvals[i]+xvals[i+1]);
		EXPECT_NEAR( interp1->f(x), ymids[i], tolerance*std::fabs(ymids[i]) );
	}
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ComplexInterpolatedValuesAreCorrect) {
	for (size_t i = 0; i < 4; i++) {
		double x = 0.5*(xvals[i]+xvals[i+1]);
		EXPECT_NEAR( c_interp1->cf(x).real(), ymids[i], tolerance*std::fabs(ymids[i]) );
		EXPECT_NEAR( c_interp1->cf(x).imag(), imids[i], tolerance*std::fabs(imids[i]) );
	}
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,RealInterpolatedFirstDerivativesAreCorrect) {
	for (size_t i = 0; i < 4; i++) {
		double x = 0.5*(xvals[i]+xvals[i+1]);
		EXPECT_NEAR( interp1->df(x), derivs1[i], tolerance*std::fabs(derivs1[i]) );
	}
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ComplexInterpolatedFirstDerivativesAreCorrect) {
	for (size_t i = 0; i < 4; i++) {
		double x = 0.5*(xvals[i]+xvals[i+1]);
		EXPECT_NEAR( c_interp1->cdf(x).real(), derivs1[i], tolerance*std::fabs(derivs1[i]) );
		EXPECT_NEAR( c_interp1->cdf(x).imag(), iderivs1[i], tolerance*std::fabs(iderivs1[i]) );
	}
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ThrowsDomainErrorIfNotReady) {
	interp1->set( 5, xvals, cvals );
	EXPECT_FALSE( interp1->is_ready() );
	EXPECT_THROW( {double d = interp1->df(2,1.5);}, std::domain_error );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ReturnsCorrectMaxDerivativeOrder) {
	EXPECT_EQ(interp1->max_derivative(),3);
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,IsExtrapolatingReturnsTrueByDefault) {
	EXPECT_TRUE( interp1->is_extrapolating() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,IsExtrapolatingReturnsFalseIfSet) {
	interp1->is_extrapolating( false );
	EXPECT_FALSE( interp1->is_extrapolating() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ExtrapolationReturnsCorrectRealValueOnLowEnd) {
	EXPECT_DOUBLE_EQ( interp1->f(-1.0), interp1->f(1.0) - 2.0*interp1->df(1.0) );
	EXPECT_DOUBLE_EQ( interp1->f(-3.0), interp1->f(1.0) - 4.0*interp1->df(1.0) );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ExtrapolationReturnsCorrectRealValueOnHighEnd) {
	EXPECT_DOUBLE_EQ( interp1->f(7.0), interp1->f(5.0) + 2.0*interp1->df(5.0) );
	EXPECT_DOUBLE_EQ( interp1->f(9.0), interp1->f(5.0) + 4.0*interp1->df(5.0) );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ExtrapolationReturnsCorrectComplexValueOnLowEnd) {
	EXPECT_DOUBLE_EQ( c_interp1->cf(-1.0).real(),
			c_interp1->cf(1.0).real() - 2.0*c_interp1->cdf(1.0).real() );
	EXPECT_DOUBLE_EQ( c_interp1->cf(-1.0).imag(),
			c_interp1->cf(1.0).imag() - 2.0*c_interp1->cdf(1.0).imag() );
	EXPECT_DOUBLE_EQ( c_interp1->cf(-4.0).real(),
			c_interp1->cf(1.0).real() - 5.0*c_interp1->cdf(1.0).real() );
	EXPECT_DOUBLE_EQ( c_interp1->cf(-4.0).imag(),
			c_interp1->cf(1.0).imag() - 5.0*c_interp1->cdf(1.0).imag() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ExtrapolationReturnsCorrectComplexValueOnHighEnd) {
	EXPECT_DOUBLE_EQ( c_interp1->cf(7.0).real(),
			c_interp1->cf(5.0).real() + 2.0*c_interp1->cdf(5.0).real() );
	EXPECT_DOUBLE_EQ( c_interp1->cf(7.0).imag(),
			c_interp1->cf(5.0).imag() + 2.0*c_interp1->cdf(5.0).imag() );
	EXPECT_DOUBLE_EQ( c_interp1->cf(9.0).real(),
			c_interp1->cf(5.0).real() + 4.0*c_interp1->cdf(5.0).real() );
	EXPECT_DOUBLE_EQ( c_interp1->cf(9.0).imag(),
			c_interp1->cf(5.0).imag() + 4.0*c_interp1->cdf(5.0).imag() );
}

TEST_F(LANLNaturalCubicSplineInterpolator1DTest,ExtrapolationThrowsOutOfRangeIfExtrapolationDisabled) {
	interp1->is_extrapolating( false );
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
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
