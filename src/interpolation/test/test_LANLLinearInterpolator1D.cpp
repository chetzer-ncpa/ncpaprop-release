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
class LANLLinearInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		interp1 = Interpolator1D::build( interpolator1d_t::LANL_1D_LINEAR );
		interp1->set( 5, xvals, yvals )->ready();

		c_interp1 = Interpolator1D::build( interpolator1d_t::LANL_1D_LINEAR );
		c_interp1->set( 5, xvals, cvals )->ready();

		pInterp1 = static_cast<LANLLinearInterpolator1D*>(interp1);
		pCInterp1 = static_cast<LANLLinearInterpolator1D*>(c_interp1);
	}

	// If there's any cleanup other than normal destructors
	void TearDown() override {
		delete interp1;
		delete c_interp1;
	}

	// class members here
	Interpolator1D *interp1, *c_interp1;
	LANLLinearInterpolator1D *pInterp1, *pCInterp1;
	double xvals[5] = { 1, 2, 3, 4, 5 };
	double yvals[5] = { 2, 4, 8, 16, 32 };
	double zvals[5] = { 0, 0, 0, 0, 0 };
	double ivals[5] = { 2, 1, 0, -2, -4 };
	complex<double> cvals[5] = {
			complex<double>( 2, 2 ),
			complex<double>( 4, 1 ),
			complex<double>( 8, 0 ),
			complex<double>( 16, -2 ),
			complex<double>( 32, -4 )
	};
};


// Tests that use the fixture
TEST_F(LANLLinearInterpolator1DTest,CopyConstructorCreatesCopy) {
	LANLLinearInterpolator1D interp2( *pInterp1 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->x_vals, pInterp1->get_real_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->f_vals, pInterp1->get_real_spline()->f_vals );
	EXPECT_EQ( interp2.is_ready(), pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,SwapWorksProperly) {
	LANLLinearInterpolator1D interp2( *pInterp1 ), c_interp2(*pCInterp1);
	swap( interp2, c_interp2 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_real_spline()->x_vals, pInterp1->get_real_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_real_spline()->f_vals, pInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->x_vals, pCInterp1->get_real_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_real_spline()->f_vals, pCInterp1->get_real_spline()->f_vals );

	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_imag_spline()->x_vals, pCInterp1->get_imag_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.get_imag_spline()->f_vals, pCInterp1->get_imag_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_imag_spline()->x_vals, pInterp1->get_imag_spline()->x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.get_imag_spline()->f_vals, pInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,RealSetOverwritesPriorValues) {
	c_interp1->set( 5, xvals, xvals )->ready();
	EXPECT_DOUBLE_ARRAY_EQ( 5, xvals, pCInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pCInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	interp1->set( 5, xvals, cvals )->ready();
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pCInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pCInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,RealSetMakesUnready) {
	c_interp1->set( 5, xvals, xvals );
	EXPECT_FALSE( c_interp1->is_ready() );
	c_interp1->ready();
	EXPECT_DOUBLE_ARRAY_EQ( 5, xvals, pCInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pCInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexSetMakesUnready) {
	interp1->set( 5, xvals, cvals );
	EXPECT_FALSE( interp1->is_ready() );
	interp1->ready();
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pCInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pCInterp1->get_imag_spline()->f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,ThrowsDomainErrorIfNotReady) {
	interp1->set( 5, xvals, cvals );
	EXPECT_FALSE( interp1->is_ready() );
	EXPECT_THROW( {double d = interp1->df(2,1.5);}, std::domain_error );
}

TEST_F(LANLLinearInterpolator1DTest,SetXOverwritesValuesAndMakesUnready) {
	pInterp1->set_x( 5, yvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pInterp1->get_real_spline()->x_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,RealSetYOverwritesValuesAndMakesUnready) {
	pInterp1->set_y( 5, zvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pInterp1->get_real_spline()->f_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexSetYOverwritesValuesAndMakesUnready) {
	pInterp1->set_y( 5, cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pInterp1->get_real_spline()->f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pInterp1->get_imag_spline()->f_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,InitClearsSplines) {
	interp1->init();
	EXPECT_EQ( pInterp1->get_real_spline()->length,0);
	EXPECT_EQ( pInterp1->get_imag_spline()->length,0);
}

TEST_F(LANLLinearInterpolator1DTest,InitSetsPointersToNull) {
	interp1->init();
	EXPECT_EQ( pInterp1->get_real_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_real_spline()->f_vals, nullptr);
	EXPECT_EQ( pInterp1->get_imag_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_imag_spline()->f_vals, nullptr);
}

TEST_F(LANLLinearInterpolator1DTest,AllocateSetsPointersToValid) {
	interp1->init();
	interp1->allocate( 15 );
	EXPECT_NE( pInterp1->get_real_spline()->x_vals, nullptr );
	EXPECT_NE( pInterp1->get_real_spline()->f_vals, nullptr);
	EXPECT_NE( pInterp1->get_imag_spline()->x_vals, nullptr );
	EXPECT_NE( pInterp1->get_imag_spline()->f_vals, nullptr);
}

TEST_F(LANLLinearInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
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

TEST_F(LANLLinearInterpolator1DTest,FreeSetsPointersToNull) {
	interp1->free();
	EXPECT_EQ( pInterp1->get_real_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_real_spline()->f_vals, nullptr);
	EXPECT_EQ( pInterp1->get_imag_spline()->x_vals, nullptr );
	EXPECT_EQ( pInterp1->get_imag_spline()->f_vals, nullptr);
}

TEST_F(LANLLinearInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( interp1->is_ready() );
	pInterp1->free();
	EXPECT_FALSE( interp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,IdentifierIsCorrect) {
	EXPECT_EQ( interp1->identifier(), "LANL 1-D Linear Interpolator");
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedValuesAreCorrect) {
	EXPECT_NEAR( interp1->f( 1.5 ), 3.0, 1e-5 );
	EXPECT_NEAR( interp1->f( 2.5 ), 6.0, 1e-5 );
	EXPECT_NEAR( interp1->f( 3.5 ), 12.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolatedValuesAreCorrect) {
	EXPECT_NEAR( c_interp1->cf( 1.5 ).real(), 3.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 2.5 ).real(), 6.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 3.5 ).real(), 12.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 1.5 ).imag(), 1.5, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 2.5 ).imag(), 0.5, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 3.5 ).imag(), -1.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedFirstDerivativesAreCorrect) {
	EXPECT_NEAR( interp1->df( 1.5 ), 2.0, 1e-5 );
	EXPECT_NEAR( interp1->df( 2.5 ), 4.0, 1e-5 );
	EXPECT_NEAR( interp1->df( 3.5 ), 8.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolatedFirstDerivativesAreCorrect) {
	EXPECT_NEAR( c_interp1->cdf( 1.5 ).real(), 2.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 2.5 ).real(), 4.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 3.5 ).real(), 8.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 1.5 ).imag(), -1.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 2.5 ).imag(), -1.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 3.5 ).imag(), -2.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedValuesEqualDerivativesCallWithN0) {
	EXPECT_DOUBLE_EQ( interp1->df( 0, 1.5 ), interp1->f( 1.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 0, 2.5 ), interp1->f( 2.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 0, 3.5 ), interp1->f( 3.5 ) );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolatedValuesEqualDerivativesCallWithN0) {
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 1.5 ).real(), interp1->cf( 1.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 2.5 ).real(), interp1->cf( 2.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 3.5 ).real(), interp1->cf( 3.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 1.5 ).imag(), interp1->cf( 1.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 2.5 ).imag(), interp1->cf( 2.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 3.5 ).imag(), interp1->cf( 3.5 ).imag() );
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedFirstDerivativesEqualCallWithN1) {
	EXPECT_DOUBLE_EQ( interp1->df( 1, 1.5 ), interp1->df( 1.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 1, 2.5 ), interp1->df( 2.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 1, 3.5 ), interp1->df( 3.5 ) );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolatedFirstDerivativesEqualCallWithN1) {
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 1.5 ).real(), interp1->cdf( 1.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 2.5 ).real(), interp1->cdf( 2.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 3.5 ).real(), interp1->cdf( 3.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 1.5 ).imag(), interp1->cdf( 1.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 2.5 ).imag(), interp1->cdf( 2.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 3.5 ).imag(), interp1->cdf( 3.5 ).imag() );
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedSecondDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( interp1->df( 2, 1.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 2, 2.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 2, 3.5 ), 0.0 );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolateSecondDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 1.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 2.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 3.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 1.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 2.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 3.5 ).imag(), 0.0 );
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedThirdDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( interp1->df( 3, 1.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 3, 2.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 3, 3.5 ), 0.0 );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolateThirdDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 1.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 2.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 3.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 1.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 2.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 3.5 ).imag(), 0.0 );
}

TEST_F(LANLLinearInterpolator1DTest,IsExtrapolatingReturnsTrueByDefault) {
	EXPECT_TRUE( interp1->is_extrapolating() );
}

TEST_F(LANLLinearInterpolator1DTest,IsExtrapolatingReturnsFalseIfSet) {
	interp1->is_extrapolating( false );
	EXPECT_FALSE( interp1->is_extrapolating() );
}

TEST_F(LANLLinearInterpolator1DTest,ExtrapolationThrowsOutOfRangeIfExtrapolationDisabled) {
	interp1->is_extrapolating( false );
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
}


TEST_F(LANLLinearInterpolator1DTest,ReturnsCorrectMaxDerivativeOrder) {
	EXPECT_EQ(interp1->max_derivative(),3);
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
