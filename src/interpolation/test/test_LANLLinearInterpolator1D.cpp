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
		interp1->set( 5, xvals, yvals );

		c_interp1 = Interpolator1D::build( interpolator1d_t::LANL_1D_LINEAR );
		c_interp1->set( 5, xvals, cvals );

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
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.real_spline_.x_vals, pInterp1->real_spline_.x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.real_spline_.f_vals, pInterp1->real_spline_.f_vals );
	EXPECT_EQ( interp2.ready_, pInterp1->ready_ );
}

TEST_F(LANLLinearInterpolator1DTest,SwapWorksProperly) {
	LANLLinearInterpolator1D interp2( *pInterp1 ), c_interp2(*pCInterp1);
	swap( interp2, c_interp2 );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.real_spline_.x_vals, pInterp1->real_spline_.x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.real_spline_.f_vals, pInterp1->real_spline_.f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.real_spline_.x_vals, pCInterp1->real_spline_.x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.real_spline_.f_vals, pCInterp1->real_spline_.f_vals );

	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.imag_spline_.x_vals, pCInterp1->imag_spline_.x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, interp2.imag_spline_.f_vals, pCInterp1->imag_spline_.f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.imag_spline_.x_vals, pInterp1->imag_spline_.x_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, c_interp2.imag_spline_.f_vals, pInterp1->imag_spline_.f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,RealSetOverwritesPriorValues) {
	c_interp1->set( 5, xvals, xvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, xvals, pCInterp1->real_spline_.f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pCInterp1->imag_spline_.f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexSetOverwritesPriorValues) {
	interp1->set( 5, xvals, cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pCInterp1->real_spline_.f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pCInterp1->imag_spline_.f_vals );
}

TEST_F(LANLLinearInterpolator1DTest,SetXOverwritesValuesAndMakesUnready) {
	pInterp1->set_x( 5, yvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pInterp1->real_spline_.x_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,RealSetYOverwritesValuesAndMakesUnready) {
	pInterp1->set_y( 5, zvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, zvals, pInterp1->real_spline_.f_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexSetYOverwritesValuesAndMakesUnready) {
	pInterp1->set_y( 5, cvals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, yvals, pInterp1->real_spline_.f_vals );
	EXPECT_DOUBLE_ARRAY_EQ( 5, ivals, pInterp1->imag_spline_.f_vals );
	EXPECT_FALSE( pInterp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,InitClearsSplines) {
	interp1->init();
	EXPECT_EQ( pInterp1->real_spline_.length,0);
	EXPECT_EQ( pInterp1->imag_spline_.length,0);
}

TEST_F(LANLLinearInterpolator1DTest,InitSetsPointersToNull) {
	interp1->init();
	EXPECT_EQ( pInterp1->real_spline_.x_vals, nullptr );
	EXPECT_EQ( pInterp1->real_spline_.f_vals, nullptr);
	EXPECT_EQ( pInterp1->imag_spline_.x_vals, nullptr );
	EXPECT_EQ( pInterp1->imag_spline_.f_vals, nullptr);
}

TEST_F(LANLLinearInterpolator1DTest,AllocateSetsPointersToValid) {
	interp1->init();
	interp1->allocate( 15 );
	EXPECT_NE( pInterp1->real_spline_.x_vals, nullptr );
	EXPECT_NE( pInterp1->real_spline_.f_vals, nullptr);
	EXPECT_NE( pInterp1->imag_spline_.x_vals, nullptr );
	EXPECT_NE( pInterp1->imag_spline_.f_vals, nullptr);
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
	EXPECT_EQ( pInterp1->real_spline_.x_vals, nullptr );
	EXPECT_EQ( pInterp1->real_spline_.f_vals, nullptr);
	EXPECT_EQ( pInterp1->imag_spline_.x_vals, nullptr );
	EXPECT_EQ( pInterp1->imag_spline_.f_vals, nullptr);
}

TEST_F(LANLLinearInterpolator1DTest,IsReadyReturnsProperValue) {
	EXPECT_TRUE( interp1->is_ready() );
	pInterp1->ready_ = false;
	EXPECT_FALSE( interp1->is_ready() );
}

TEST_F(LANLLinearInterpolator1DTest,IdentifierIsCorrect) {
	EXPECT_EQ( interp1->identifier(), "LANL 1-D Linear Interpolator");
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedValuesAreCorrect) {
	EXPECT_NEAR( interp1->eval_f( 1.5 ), 3.0, 1e-5 );
	EXPECT_NEAR( interp1->eval_f( 2.5 ), 6.0, 1e-5 );
	EXPECT_NEAR( interp1->eval_f( 3.5 ), 12.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolatedValuesAreCorrect) {
	EXPECT_NEAR( c_interp1->eval_cf( 1.5 ).real(), 3.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cf( 2.5 ).real(), 6.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cf( 3.5 ).real(), 12.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cf( 1.5 ).imag(), 1.5, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cf( 2.5 ).imag(), 0.5, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cf( 3.5 ).imag(), -1.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,RealInterpolatedFirstDerivativesAreCorrect) {
	EXPECT_NEAR( interp1->eval_df( 1.5 ), 2.0, 1e-5 );
	EXPECT_NEAR( interp1->eval_df( 2.5 ), 4.0, 1e-5 );
	EXPECT_NEAR( interp1->eval_df( 3.5 ), 8.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,ComplexInterpolatedFirstDerivativesAreCorrect) {
	EXPECT_NEAR( c_interp1->eval_cdf( 1.5 ).real(), 2.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cdf( 2.5 ).real(), 4.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cdf( 3.5 ).real(), 8.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cdf( 1.5 ).imag(), -1.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cdf( 2.5 ).imag(), -1.0, 1e-5 );
	EXPECT_NEAR( c_interp1->eval_cdf( 3.5 ).imag(), -2.0, 1e-5 );
}

TEST_F(LANLLinearInterpolator1DTest,ThrowsDomainErrorForSecondDerivative) {
	EXPECT_THROW( {double d = interp1->eval_df(2,1.5);}, std::domain_error );
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
