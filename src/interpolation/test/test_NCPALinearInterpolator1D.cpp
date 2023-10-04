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
class NCPALinearInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_LINEAR );
		interp1->set( 5, xvals, yvals )->ready();

		c_interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_LINEAR );
		c_interp1->set( 5, xvals, cvals )->ready();

		pInterp1 = static_cast<NCPALinearInterpolator1D*>(interp1);
		pCInterp1 = static_cast<NCPALinearInterpolator1D*>(c_interp1);
	}

	// If there's any cleanup other than normal destructors
	void TearDown() override {
		delete interp1;
		delete c_interp1;
	}

	// class members here
	Interpolator1D *interp1, *c_interp1;
	NCPALinearInterpolator1D *pInterp1, *pCInterp1;
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
TEST_F(NCPALinearInterpolator1DTest,CopyConstructorCreatesCopy) {
	NCPALinearInterpolator1D interp2( *pInterp1 );
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( interp2.f(xvals[i]), pInterp1->f(xvals[i]));
	}
	EXPECT_EQ( interp2.is_ready(), pInterp1->is_ready() );
}

TEST_F(NCPALinearInterpolator1DTest,SwapWorksProperly) {
	NCPALinearInterpolator1D interp2( *pInterp1 ), c_interp2(*pCInterp1);
	swap( interp2, c_interp2 );
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( c_interp2.f(xvals[i]), pInterp1->f(xvals[i]));
		EXPECT_DOUBLE_EQ( interp2.f(xvals[i]), pCInterp1->f(xvals[i]));
	}
}

TEST_F(NCPALinearInterpolator1DTest,RealSetOverwritesPriorValues) {
	c_interp1->set( 5, xvals, xvals )->ready();
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( xvals[i], c_interp1->f(xvals[i]));
	}
}

TEST_F(NCPALinearInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	interp1->set( 5, xvals, cvals )->ready();
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( yvals[i], pCInterp1->cf(xvals[i]).real());
		EXPECT_DOUBLE_EQ( ivals[i], pCInterp1->cf(xvals[i]).imag());
	}
}

TEST_F(NCPALinearInterpolator1DTest,ThrowsDomainErrorIfNotReady) {
	interp1->init();
	EXPECT_THROW( {double d = interp1->df(2,1.5);}, std::domain_error );
}

TEST_F(NCPALinearInterpolator1DTest,InitClearsSplines) {
	interp1->init();
	EXPECT_FALSE( interp1->is_ready() );
}


TEST_F(NCPALinearInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
	interp1->init();
	EXPECT_FALSE( interp1->is_ready() );
	interp1->allocate(5);
	EXPECT_FALSE( interp1->is_ready() );
	pInterp1->set( 5, yvals, zvals );
	EXPECT_TRUE( interp1->is_ready() );
	interp1->ready();
	EXPECT_TRUE( interp1->is_ready() );
}


TEST_F(NCPALinearInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( interp1->is_ready() );
	pInterp1->free();
	EXPECT_FALSE( interp1->is_ready() );
}

TEST_F(NCPALinearInterpolator1DTest,IdentifierIsCorrect) {
	EXPECT_EQ( interp1->identifier(), Interpolator1D::as_string( interp1->type() ));
}

TEST_F(NCPALinearInterpolator1DTest,RealInterpolatedValuesAreCorrect) {
	EXPECT_NEAR( interp1->f( 1.5 ), 3.0, 1e-5 );
	EXPECT_NEAR( interp1->f( 2.5 ), 6.0, 1e-5 );
	EXPECT_NEAR( interp1->f( 3.5 ), 12.0, 1e-5 );
}

TEST_F(NCPALinearInterpolator1DTest,ComplexInterpolatedValuesAreCorrect) {
	EXPECT_NEAR( c_interp1->cf( 1.5 ).real(), 3.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 2.5 ).real(), 6.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 3.5 ).real(), 12.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 1.5 ).imag(), 1.5, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 2.5 ).imag(), 0.5, 1e-5 );
	EXPECT_NEAR( c_interp1->cf( 3.5 ).imag(), -1.0, 1e-5 );
}

TEST_F(NCPALinearInterpolator1DTest,RealInterpolatedFirstDerivativesAreCorrect) {
	EXPECT_NEAR( interp1->df( 1.5 ), 2.0, 1e-5 );
	EXPECT_NEAR( interp1->df( 2.5 ), 4.0, 1e-5 );
	EXPECT_NEAR( interp1->df( 3.5 ), 8.0, 1e-5 );
}

TEST_F(NCPALinearInterpolator1DTest,ComplexInterpolatedFirstDerivativesAreCorrect) {
	EXPECT_NEAR( c_interp1->cdf( 1.5 ).real(), 2.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 2.5 ).real(), 4.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 3.5 ).real(), 8.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 1.5 ).imag(), -1.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 2.5 ).imag(), -1.0, 1e-5 );
	EXPECT_NEAR( c_interp1->cdf( 3.5 ).imag(), -2.0, 1e-5 );
}

TEST_F(NCPALinearInterpolator1DTest,RealInterpolatedValuesEqualDerivativesCallWithN0) {
	EXPECT_DOUBLE_EQ( interp1->df( 0, 1.5 ), interp1->f( 1.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 0, 2.5 ), interp1->f( 2.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 0, 3.5 ), interp1->f( 3.5 ) );
}

TEST_F(NCPALinearInterpolator1DTest,ComplexInterpolatedValuesEqualDerivativesCallWithN0) {
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 1.5 ).real(), interp1->cf( 1.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 2.5 ).real(), interp1->cf( 2.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 3.5 ).real(), interp1->cf( 3.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 1.5 ).imag(), interp1->cf( 1.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 2.5 ).imag(), interp1->cf( 2.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 3.5 ).imag(), interp1->cf( 3.5 ).imag() );
}

TEST_F(NCPALinearInterpolator1DTest,RealInterpolatedFirstDerivativesEqualCallWithN1) {
	EXPECT_DOUBLE_EQ( interp1->df( 1, 1.5 ), interp1->df( 1.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 1, 2.5 ), interp1->df( 2.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 1, 3.5 ), interp1->df( 3.5 ) );
}

TEST_F(NCPALinearInterpolator1DTest,ComplexInterpolatedFirstDerivativesEqualCallWithN1) {
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 1.5 ).real(), interp1->cdf( 1.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 2.5 ).real(), interp1->cdf( 2.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 3.5 ).real(), interp1->cdf( 3.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 1.5 ).imag(), interp1->cdf( 1.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 2.5 ).imag(), interp1->cdf( 2.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 3.5 ).imag(), interp1->cdf( 3.5 ).imag() );
}

TEST_F(NCPALinearInterpolator1DTest,RealInterpolatedSecondDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( interp1->df( 2, 1.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 2, 2.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 2, 3.5 ), 0.0 );
}

TEST_F(NCPALinearInterpolator1DTest,ComplexInterpolateSecondDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 1.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 2.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 3.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 1.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 2.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 3.5 ).imag(), 0.0 );
}

TEST_F(NCPALinearInterpolator1DTest,RealInterpolatedThirdDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( interp1->df( 3, 1.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 3, 2.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 3, 3.5 ), 0.0 );
}

TEST_F(NCPALinearInterpolator1DTest,ComplexInterpolateThirdDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 1.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 2.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 3.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 1.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 2.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 3.5 ).imag(), 0.0 );
}

TEST_F(NCPALinearInterpolator1DTest,IsExtrapolatingReturnsTrueByDefault) {
	EXPECT_TRUE( interp1->is_extrapolating() );
}

TEST_F(NCPALinearInterpolator1DTest,IsExtrapolatingReturnsFalseIfSet) {
	interp1->is_extrapolating( false );
	EXPECT_FALSE( interp1->is_extrapolating() );
}

TEST_F(NCPALinearInterpolator1DTest,ExtrapolationThrowsOutOfRangeIfExtrapolationDisabled) {
	interp1->is_extrapolating( false );
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
}

TEST_F(NCPALinearInterpolator1DTest,ExtrapolationReturnsCorrectRealValueOnLowEnd) {
	EXPECT_DOUBLE_EQ( interp1->f(-1.0), -2.0 );
	EXPECT_DOUBLE_EQ( interp1->f(-3.0), -6.0 );
}

TEST_F(NCPALinearInterpolator1DTest,ExtrapolationReturnsCorrectRealValueOnHighEnd) {
	EXPECT_DOUBLE_EQ( interp1->f(7.0), 64.0 );
	EXPECT_DOUBLE_EQ( interp1->f(9.0), 96.0 );
}

TEST_F(NCPALinearInterpolator1DTest,ExtrapolationReturnsCorrectComplexValueOnLowEnd) {
	EXPECT_DOUBLE_EQ( c_interp1->cf(-1.0).real(), -2.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cf(-1.0).imag(), 4.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cf(-4.0).real(), -8.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cf(-4.0).imag(), 7.0 );
}

TEST_F(NCPALinearInterpolator1DTest,ExtrapolationReturnsCorrectComplexValueOnHighEnd) {
	EXPECT_DOUBLE_EQ( c_interp1->cf(7.0).real(), 64.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cf(7.0).imag(), -8.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cf(9.0).real(), 96.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cf(9.0).imag(), -12.0 );
}


TEST_F(NCPALinearInterpolator1DTest,ReturnsCorrectMaxDerivativeOrder) {
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
