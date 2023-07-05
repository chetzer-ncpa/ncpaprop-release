#include "util.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace std;
using namespace NCPA;
using namespace testing;


// Test close_enough function
TEST(utilTest, CloseEnoughHandlesPositiveInputs) {
	EXPECT_TRUE( within( 2.5366672, 2.5366685, 4) );
	EXPECT_FALSE( within( 2.5366672, 2.5366685, 6) );
}

TEST(utilTest, CloseEnoughHandlesNegativeInputs) {
	EXPECT_TRUE( within( -2.5366672, -2.5366685, 4) );
    EXPECT_FALSE( within( -2.5366672, -2.5366685, 6) );

}

TEST(utilTest, CircShiftShiftsPositively) {
	size_t invec3[3] = { 1, 2, 3 }, outvec3[3];
	circshift<size_t>(invec3,3,2,outvec3);
	ASSERT_THAT(outvec3,ElementsAre(3, 1, 2));
	size_t invec5[5] = { 1, 2, 3, 4, 5 }, outvec5[5];
	circshift<size_t>(invec5,5,2,outvec5);
	ASSERT_THAT(outvec5,ElementsAre(3, 4, 5, 1, 2));
}

TEST(utilTest, CircShiftShiftsNegatively) {
	size_t invec3[3] = { 1, 2, 3 }, outvec3[3];
	circshift<size_t>(invec3,3,-1,outvec3);
	ASSERT_THAT(outvec3,ElementsAre(3, 1, 2));
	size_t invec5[5] = { 1, 2, 3, 4, 5 }, outvec5[5];
	circshift<size_t>(invec5,5,-2,outvec5);
	ASSERT_THAT(outvec5,ElementsAre(4, 5, 1, 2, 3));
}

TEST(utilTest,FillMatrixFillsRowVector) {
	double **mat = NCPA::allocate_matrix<double>(3,1);
	EXPECT_DOUBLE_EQ( mat[0][0], 0.0 );
	NCPA::fill_matrix<double>(mat,3,1,6.0);
	EXPECT_DOUBLE_EQ( mat[2][0], 6.0 );
}

TEST(utilTest,FindClosestIndexFindsInternalMatch) {
	double testvec[ 6 ] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0};
	double testval = 5.0;
	ASSERT_EQ( NCPA::find_closest_index(testvec, 6, testval), 2 );
}

TEST(utilTest,FindClosestIndexFindsEndMatches) {
	double testvec[ 6 ] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0};
	double lowtestval = -1.0, hightestval = 19.5;
	ASSERT_EQ( NCPA::find_closest_index(testvec, 6, lowtestval), 0 );
	ASSERT_EQ( NCPA::find_closest_index(testvec, 6, hightestval), 5 );
}

TEST(utilTest,Cart2PolConvertsCorrectly) {
	double rho, theta, x, y;

	// test 1: [0,4]
	x = 0; y = 4;
	NCPA::cart2pol(x,y,rho,theta);
	ASSERT_DOUBLE_EQ( rho, 4.0 );
	ASSERT_DOUBLE_EQ( theta, M_PI / 2.0 );

	// test 2: [4,0]
	x = 4; y = 0;
	NCPA::cart2pol(x,y,rho,theta);
	ASSERT_DOUBLE_EQ( rho, 4.0 );
	ASSERT_DOUBLE_EQ( theta, 0.0 );

	// test 3: 4, -4
	x = 4; y = -4;
	NCPA::cart2pol(x,y,rho,theta);
	ASSERT_DOUBLE_EQ( rho, std::sqrt( 32.0 ) );
	ASSERT_DOUBLE_EQ( theta, -M_PI / 4.0 );
}

TEST(utilTest,EvalPolyReturnsConstantTermWithZeroInput) {
	double c[ 4 ] = { 5.0, 2.0, 3.4, 1.2 };
	double x = 0.0; double expected = c[0];
	ASSERT_DOUBLE_EQ( NCPA::evalpoly( 4, c, x ), expected );
}

TEST(utilTest,EvalPolyCorrectWithPositiveInput) {
	double c[ 4 ] = { 5.0, -2.0, 3.4, 1.2 };
	double x = M_PI; double expected = c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x;
	ASSERT_DOUBLE_EQ( NCPA::evalpoly( 4, c, x ), expected );
}

TEST(utilTest,EvalPolyCorrectWithNegativeInput) {
	double c[ 4 ] = { 5.0, -2.0, 3.4, -1.2 };
	double x = -M_PI / 4.0; double expected = c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x;
	ASSERT_DOUBLE_EQ( NCPA::evalpoly( 4, c, x ), expected );
}

TEST(utilTest,Double2ComplexWithRealOnlyReturnsCorrectValues) {
	double in[ 5 ] = { 1, -1, 2, -2, 5 };
	complex<double> out[5];

	double2complex( 5, in, out );
	for (size_t i = 0; i < 5; i++) {
		ASSERT_DOUBLE_EQ( out[i].real(), in[i] );
		ASSERT_DOUBLE_EQ( out[i].imag(), 0.0 );
	}
}

TEST(utilTest,Double2ComplexReturnsCorrectValues) {
	double real[ 5 ] = { 1, -1, 2, -2, 5 };
	double imag[ 5 ] = { 2, 4, 6, 8, 10 };
	complex<double> out[5];

	double2complex( 5, real, imag, out );
	for (size_t i = 0; i < 5; i++) {
		ASSERT_DOUBLE_EQ( out[i].real(), real[i] );
		ASSERT_DOUBLE_EQ( out[i].imag(), imag[i] );
	}
}

TEST(utilTest,Complex2DoubleReturnsCorrectValues) {
	double real[ 5 ] = { 1, -1, 2, -2, 5 };
	double imag[ 5 ] = { 2, 4, 6, 8, 10 };
	complex<double> in[5];
	double out_r[5], out_i[5];
	for (size_t i = 0; i < 5; i++) {
		in[i] = complex<double>( real[i], imag[i] );
	}

	complex2double( 5, in, out_r, out_i );
	for (size_t i = 0; i < 5; i++) {
		ASSERT_DOUBLE_EQ( out_r[i], real[i] );
		ASSERT_DOUBLE_EQ( out_i[i], imag[i] );
	}
}


