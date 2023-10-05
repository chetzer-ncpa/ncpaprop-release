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
class Interpolator1DTest : public ::testing::Test {
protected:
	void SetUp() override {
		xvector.resize(5);
		yvector.resize(5);
		cvector.resize(5);
		std::copy( xvals, xvals+5, xvector.begin() );
		std::copy( yvals, yvals+5, yvector.begin() );
		std::copy( cvals, cvals+5, cvector.begin() );
		v_interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_LINEAR );
		v_interp1->set( xvector, yvector )->ready();

		vc_interp1 = Interpolator1D::build( interpolator1d_t::NCPA_1D_LINEAR );
		vc_interp1->set( xvector, cvector )->ready();
	}

	void TearDown() override {
		delete v_interp1;
		delete vc_interp1;
	}

	Interpolator1D *iPtr;
	std::vector<interpolator1d_t> types = {
			interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR,
			interpolator1d_t::NCPA_1D_LINEAR,
			interpolator1d_t::LANL_1D_LINEAR,
			interpolator1d_t::LANL_1D_NATURAL_CUBIC,
			interpolator1d_t::GSL_1D_LINEAR,
			interpolator1d_t::GSL_1D_POLYNOMIAL,
			interpolator1d_t::GSL_1D_CUBIC,
			interpolator1d_t::GSL_1D_CUBIC_PERIODIC,
			interpolator1d_t::GSL_1D_AKIMA,
			interpolator1d_t::GSL_1D_AKIMA_PERIODIC,
			interpolator1d_t::GSL_1D_STEFFEN
	};
	Interpolator1D *v_interp1, *vc_interp1;
	double xvals[5] = { 1, 2, 3, 4, 5 };
	double yvals[5] = { 2, 4, 8, 16, 32 };
	std::complex<double> cvals[5] = {
			complex<double>( 2, 2 ),
			complex<double>( 4, 1 ),
			complex<double>( 8, 0 ),
			complex<double>( 16, -2 ),
			complex<double>( 32, -4 )
	};
	std::vector<double> yvector, xvector;
	std::vector<std::complex<double>> cvector;

};


TEST_F(Interpolator1DTest,CanBuildReturnsTrueForNearestNeighbor) {
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::NCPA_1D_NEAREST_NEIGHBOR ) );
}

TEST_F(Interpolator1DTest,CanBuildReturnsTrueForNCPALinear) {
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::NCPA_1D_LINEAR ) );
}

TEST_F(Interpolator1DTest,CanBuildReturnsCorrectBooleanForLANL) {
#ifdef HAVE_LANL_INTERPOLATION_LIBRARY
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::LANL_1D_LINEAR ) );
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::LANL_1D_NATURAL_CUBIC ) );
#else
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::LANL_1D_LINEAR ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::LANL_1D_NATURAL_CUBIC ) );
#endif
}

TEST_F(Interpolator1DTest,CanBuildReturnsCorrectBooleanForGSL) {
#ifdef HAVE_GSL_INTERPOLATION_LIBRARY
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_LINEAR ) );
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_POLYNOMIAL ) );
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_CUBIC ) );
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_CUBIC_PERIODIC ) );
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_AKIMA ) );
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_AKIMA_PERIODIC ) );
#if GSL_MAJOR_VERSION > 1
	ASSERT_TRUE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_STEFFEN ) );
#else
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_STEFFEN ) );
#endif
#else
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_LINEAR ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_POLYNOMIAL ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_CUBIC ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_CUBIC_PERIODIC ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_AKIMA ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_AKIMA_PERIODIC ) );
	ASSERT_FALSE(Interpolator1D::can_build( interpolator1d_t::GSL_1D_STEFFEN ) );
#endif
}

void testInterpolatorBuild( interpolator1d_t t ) {
	Interpolator1D *iPtr;
	if (Interpolator1D::can_build( t ) ) {
		iPtr = Interpolator1D::build( t );
		ASSERT_EQ( iPtr->identifier(), Interpolator1D::as_string( iPtr->type() ) );
		delete iPtr;
	} else {
		EXPECT_THROW( {iPtr = Interpolator1D::build(t);}, std::out_of_range );
	}
}

TEST_F(Interpolator1DTest,BuildsOrThrowsAsAppropriate) {
	for (auto it = types.cbegin(); it != types.cend(); ++it) {
		testInterpolatorBuild( *it );
	}
}

TEST_F(Interpolator1DTest,RealInterpolatedValuesAreCorrectWithVectorInput) {
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( v_interp1->f( xvals[i] ), yvals[i] );
	}
}


TEST_F(Interpolator1DTest,ComplexInterpolatedValuesAreCorrectWithVectorInput) {
	for (size_t i = 0; i < 5; i++) {
		EXPECT_DOUBLE_EQ( vc_interp1->cf( xvals[i] ).real(), cvals[i].real() );
		EXPECT_DOUBLE_EQ( vc_interp1->cf( xvals[i] ).imag(), cvals[i].imag() );
	}
}


