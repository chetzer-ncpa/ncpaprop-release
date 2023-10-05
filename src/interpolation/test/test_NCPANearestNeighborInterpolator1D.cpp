#include "NCPAInterpolation.h"
#include "NCPACommon.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <complex>
#include <vector>
#include <algorithm>

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
class NCPANearestNeighborInterpolator1DTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
//		xvector.resize(5);
//		yvector.resize(5);
//		cvector.resize(5);
//		std::copy( xvals, xvals+5, xvector.begin() );
//		std::copy( yvals, yvals+5, yvector.begin() );
//		std::copy( cvals, cvals+5, cvector.begin() );

		interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR );
		interp1->set( 5, xvals, yvals )->ready();

		c_interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR );
		c_interp1->set( 5, xvals, cvals )->ready();

//		v_interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR );
//		v_interp1->set( xvector, yvector )->ready();
//
//		vc_interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR );
//		vc_interp1->set( xvector, cvector )->ready();

		pInterp1 = static_cast<NCPANearestNeighborInterpolator1D*>(interp1);
		pCInterp1 = static_cast<NCPANearestNeighborInterpolator1D*>(c_interp1);
	}

	// If there's any cleanup other than normal destructors
	void TearDown() override {
		delete interp1;
		delete c_interp1;
	}

	// class members here
	Interpolator1D *interp1, *c_interp1;
	NCPANearestNeighborInterpolator1D *pInterp1, *pCInterp1;
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
//	std::vector<double> yvector, xvector;
//	std::vector<std::complex<double>> cvector;
};


// Tests that use the fixture
TEST_F(NCPANearestNeighborInterpolator1DTest,CopyConstructorCreatesCopy) {
	NCPANearestNeighborInterpolator1D interp2( *pInterp1 );
	EXPECT_DOUBLE_EQ( interp2.get_low_interp_limit(), pInterp1->get_low_interp_limit() );
	EXPECT_DOUBLE_EQ( interp2.get_high_interp_limit(), pInterp1->get_high_interp_limit() );
	EXPECT_EQ( interp2.is_ready(), pInterp1->is_ready() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,SwapWorksProperly) {
	NCPANearestNeighborInterpolator1D interp2;
	ASSERT_FALSE( interp2.is_ready() );
	ASSERT_TRUE( pInterp1->is_ready() );
	swap( interp2, *pInterp1 );
	ASSERT_TRUE( interp2.is_ready() );
	ASSERT_FALSE( pInterp1->is_ready() );
	ASSERT_DOUBLE_EQ( interp2.f( 2.0 ), 4.0 );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,CloneWorksProperly) {
	Interpolator1D *clone = interp1->clone();
	interp1->free();
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( clone->f(xvals[i]), yvals[i] );
	}
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ThrowsDomainErrorIfNotReady) {
	interp1->init();
	EXPECT_FALSE( interp1->is_ready() );
	EXPECT_THROW( {double d = interp1->f(1.5);}, std::domain_error );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,InitClearsSplines) {
	interp1->init();
	EXPECT_FALSE( interp1->is_ready() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ReadyMakesSplinesReadyAfterSet) {
	interp1->init();
	EXPECT_FALSE( interp1->is_ready() );
	interp1->allocate(5);
	EXPECT_FALSE( interp1->is_ready() );
	pInterp1->set( 5, yvals, zvals );
	EXPECT_TRUE( interp1->is_ready() );
	interp1->ready();
	EXPECT_TRUE( interp1->is_ready() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( interp1->is_ready() );
	pInterp1->free();
	EXPECT_FALSE( interp1->is_ready() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,IdentifierIsCorrect) {
	EXPECT_EQ( interp1->identifier(), Interpolator1D::as_string( interp1->type() ));
}

TEST_F(NCPANearestNeighborInterpolator1DTest,RealInterpolatedValuesAreCorrect) {
	for (double d = 1.0; d < 4.5; d += 1.0) {
		EXPECT_DOUBLE_EQ( interp1->f( d + 0.45 ), interp1->f( d ) );
		EXPECT_DOUBLE_EQ( interp1->f( d + 0.55 ), interp1->f( d + 1.0 ) );
	}
}


TEST_F(NCPANearestNeighborInterpolator1DTest,RealExtrapolatedValuesAreCorrect) {
	EXPECT_DOUBLE_EQ( interp1->f( 0.0 ), interp1->f( 1.0 ) );
	EXPECT_DOUBLE_EQ( interp1->f( 10.0 ), interp1->f( 5.0 ) );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexInterpolatedValuesAreCorrect) {
	for (double d = 1.0; d < 4.5; d += 1.0) {
		EXPECT_DOUBLE_EQ( c_interp1->cf( d + 0.45 ).real(), c_interp1->cf( d ).real() );
		EXPECT_DOUBLE_EQ( c_interp1->cf( d + 0.45 ).imag(), c_interp1->cf( d ).imag() );
		EXPECT_DOUBLE_EQ( c_interp1->cf( d + 0.55 ).real(), c_interp1->cf( d + 1.0 ).real() );
		EXPECT_DOUBLE_EQ( c_interp1->cf( d + 0.55 ).imag(), c_interp1->cf( d + 1.0).imag() );
	}
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexExtrapolatedValuesAreCorrect) {
	EXPECT_DOUBLE_EQ( interp1->cf( 0.0 ).real(), interp1->cf( 1.0 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cf( 10.0 ).imag(), interp1->cf( 5.0 ).imag() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,RealInterpolatedFirstDerivativesAreCorrect) {
	for (double d = 1.0; d < 4.5; d += 1.0) {
		EXPECT_DOUBLE_EQ( interp1->df( d + 0.5 ), 0.0 );
	}
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexInterpolatedFirstDerivativesAreCorrect) {
	for (double d = 1.0; d < 4.5; d += 1.0) {
		EXPECT_DOUBLE_EQ( c_interp1->cdf( d + 0.5 ).real(), 0.0 );
		EXPECT_DOUBLE_EQ( c_interp1->cdf( d + 0.5 ).imag(), 0.0 );
	}
}

TEST_F(NCPANearestNeighborInterpolator1DTest,RealInterpolatedValuesEqualDerivativesCallWithN0) {
	EXPECT_DOUBLE_EQ( interp1->df( 0, 1.5 ), interp1->f( 1.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 0, 2.5 ), interp1->f( 2.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 0, 3.5 ), interp1->f( 3.5 ) );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexInterpolatedValuesEqualDerivativesCallWithN0) {
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 1.5 ).real(), interp1->cf( 1.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 2.5 ).real(), interp1->cf( 2.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 3.5 ).real(), interp1->cf( 3.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 1.5 ).imag(), interp1->cf( 1.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 2.5 ).imag(), interp1->cf( 2.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 0, 3.5 ).imag(), interp1->cf( 3.5 ).imag() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,RealInterpolatedFirstDerivativesEqualCallWithN1) {
	EXPECT_DOUBLE_EQ( interp1->df( 1, 1.5 ), interp1->df( 1.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 1, 2.5 ), interp1->df( 2.5 ) );
	EXPECT_DOUBLE_EQ( interp1->df( 1, 3.5 ), interp1->df( 3.5 ) );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexInterpolatedFirstDerivativesEqualCallWithN1) {
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 1.5 ).real(), interp1->cdf( 1.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 2.5 ).real(), interp1->cdf( 2.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 3.5 ).real(), interp1->cdf( 3.5 ).real() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 1.5 ).imag(), interp1->cdf( 1.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 2.5 ).imag(), interp1->cdf( 2.5 ).imag() );
	EXPECT_DOUBLE_EQ( interp1->cdf( 1, 3.5 ).imag(), interp1->cdf( 3.5 ).imag() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,RealInterpolatedSecondDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( interp1->df( 2, 1.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 2, 2.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 2, 3.5 ), 0.0 );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexInterpolateSecondDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 1.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 2.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 3.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 1.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 2.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 2, 3.5 ).imag(), 0.0 );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,RealInterpolatedThirdDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( interp1->df( 3, 1.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 3, 2.5 ), 0.0 );
	EXPECT_DOUBLE_EQ( interp1->df( 3, 3.5 ), 0.0 );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ComplexInterpolateThirdDerivativesAreZero) {
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 1.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 2.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 3.5 ).real(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 1.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 2.5 ).imag(), 0.0 );
	EXPECT_DOUBLE_EQ( c_interp1->cdf( 3, 3.5 ).imag(), 0.0 );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,IsExtrapolatingReturnsTrueByDefault) {
	EXPECT_TRUE( interp1->is_extrapolating() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,IsExtrapolatingSetsExtrapolatingToFalse) {
	interp1->is_extrapolating( false );
	EXPECT_FALSE( interp1->is_extrapolating() );
}

TEST_F(NCPANearestNeighborInterpolator1DTest,ExtrapolationThrowsExceptionIfExtrapolationDisabled) {
	interp1->is_extrapolating( false );
	EXPECT_THROW( {double d = interp1->f(interp1->get_low_interp_limit() - 1.0); }, std::out_of_range );
	EXPECT_THROW( {double d = interp1->f(interp1->get_high_interp_limit() + 1.0); }, std::out_of_range );
}


TEST_F(NCPANearestNeighborInterpolator1DTest,ReturnsCorrectMaxDerivativeOrder) {
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
