#include "util.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

using namespace std;
using namespace NCPA;
using namespace testing;


// Test close_enough function
TEST(closeEnoughTest, HandlesPositiveInputs) {
	EXPECT_TRUE( within( 2.5366672, 2.5366685, 4) );
	EXPECT_FALSE( within( 2.5366672, 2.5366685, 6) );
}

TEST(closeEnoughTest, HandlesNegativeInputs) {
	EXPECT_TRUE( within( -2.5366672, -2.5366685, 4) );
    EXPECT_FALSE( within( -2.5366672, -2.5366685, 6) );

}

TEST(circshiftTest, ShiftsPositively) {
	size_t invec3[3] = { 1, 2, 3 }, outvec3[3];
	circshift<size_t>(invec3,3,2,outvec3);
	ASSERT_THAT(outvec3,ElementsAre(3, 1, 2));
	size_t invec5[5] = { 1, 2, 3, 4, 5 }, outvec5[5];
	circshift<size_t>(invec5,5,2,outvec5);
	ASSERT_THAT(outvec5,ElementsAre(3, 4, 5, 1, 2));
}

TEST(circshiftTest, ShiftsNegatively) {
	size_t invec3[3] = { 1, 2, 3 }, outvec3[3];
	circshift<size_t>(invec3,3,-1,outvec3);
	ASSERT_THAT(outvec3,ElementsAre(3, 1, 2));
	size_t invec5[5] = { 1, 2, 3, 4, 5 }, outvec5[5];
	circshift<size_t>(invec5,5,-2,outvec5);
	ASSERT_THAT(outvec5,ElementsAre(4, 5, 1, 2, 3));
}

TEST(fillMatrixTest,fillRowVector) {
	double **mat = NCPA::allocate_matrix<double>(3,1);
	EXPECT_DOUBLE_EQ( mat[0][0], 0.0 );
	NCPA::fill_matrix<double>(mat,3,1,6.0);
	EXPECT_DOUBLE_EQ( mat[2][0], 6.0 );
}

TEST(findClosestIndexTest,findInternalMatch) {
	double testvec[ 6 ] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0};
	double testval = 5.0;
	ASSERT_EQ( NCPA::find_closest_index(testvec, 6, testval), 2 );
}

TEST(findClosestIndexTest,findEndMatches) {
	double testvec[ 6 ] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0};
	double lowtestval = -1.0, hightestval = 19.5;
	ASSERT_EQ( NCPA::find_closest_index(testvec, 6, lowtestval), 0 );
	ASSERT_EQ( NCPA::find_closest_index(testvec, 6, hightestval), 5 );
}

TEST(cart2polTest,testCart2PolConversions) {
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

TEST(evalPolyTest,polynomialArrayChecksZeroInput) {
	double c[ 4 ] = { 5.0, 2.0, 3.4, 1.2 };
	double x = 0.0; double expected = c[0];
	ASSERT_DOUBLE_EQ( NCPA::evalpoly( 4, c, x ), expected );
}

TEST(evalPolyTest,polynomialArrayChecksPositiveInput) {
	double c[ 4 ] = { 5.0, -2.0, 3.4, 1.2 };
	double x = M_PI; double expected = c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x;
	ASSERT_DOUBLE_EQ( NCPA::evalpoly( 4, c, x ), expected );
}

TEST(evalPolyTest,polynomialArrayChecksNegativeInput) {
	double c[ 4 ] = { 5.0, -2.0, 3.4, -1.2 };
	double x = -M_PI / 4.0; double expected = c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x;
	ASSERT_DOUBLE_EQ( NCPA::evalpoly( 4, c, x ), expected );
}
