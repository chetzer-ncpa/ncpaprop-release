#include "NCPAUnits.h"
#include "NCPAAtmosphere.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>

using namespace std;
using namespace NCPA;
using namespace testing;



// test fixture
class AtmosphericProperty1DTest : public ::testing::Test {
 protected:

	// Initializations and other setup
	void SetUp() override {
		prop2 = AtmosphericProperty1D( 10, alts, "km", temps, "C" );

		altv = VectorWithUnits(10,alts,"km");
		tempv = VectorWithUnits(10,temps,"C");
		prop3 = AtmosphericProperty1D( altv, tempv );

		prop4 = AtmosphericProperty1D( prop2 );
		prop5 = AtmosphericProperty1D( 10, alts, "km", temps_linear, "C" );
	}

	// If there's any cleanup other than normal destructors
	void TearDown() override {}

	// class members here
	AtmosphericProperty1D prop1, prop2, prop3, prop4, prop5, prop6;
	double temps[10] = { 0.0, 12.0, 8.0, 3.0, 0.0, -1.0, -1.3, -2.0, -3.0, -5.0 };
	double temps_linear[10] = { 0,1,2,3,4,5,6,7,8,9 };
	double alts[10]  = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	NCPA::VectorWithUnits altv, tempv;
};


// Tests that use the fixture
TEST_F(AtmosphericProperty1DTest,DefaultConstructorIsEmpty) {
	ASSERT_EQ(prop1.size(),0);
}

TEST_F(AtmosphericProperty1DTest,ArrayConstructorWorks) {
	ASSERT_EQ(prop2.size(),10);
	VectorWithUnits v1 = prop2.get_vector();
	VectorWithUnits a1 = prop2.get_altitude_vector();
	ASSERT_THAT(a1, ElementsAre(
				DoubleEq(0.0),
				DoubleEq(1.0),
				DoubleEq(2.0),
				DoubleEq(3.0),
				DoubleEq(4.0),
				DoubleEq(5.0),
				DoubleEq(6.0),
				DoubleEq(7.0),
				DoubleEq(8.0),
				DoubleEq(9.0)));
	ASSERT_THAT(v1, ElementsAre(
				DoubleEq(0.0),
				DoubleEq(12.0),
				DoubleEq(8.0),
				DoubleEq(3.0),
				DoubleEq(0.0),
				DoubleEq(-1.0),
				DoubleEq(-1.3),
				DoubleEq(-2.0),
				DoubleEq(-3.0),
				DoubleEq(-5.0)));
}

TEST_F(AtmosphericProperty1DTest,CopyConstructorWorks) {
	EXPECT_EQ(prop4.size(),prop2.size());
	EXPECT_EQ(prop4.get_vector(),prop2.get_vector());
	EXPECT_EQ(prop4.get_altitude_vector(),prop2.get_altitude_vector());
}

TEST_F(AtmosphericProperty1DTest,SwapWorks) {
	swap(prop2,prop6);
	EXPECT_EQ(prop6.get_altitude_vector(),prop4.get_altitude_vector());
	EXPECT_EQ(prop6.get_vector(),prop4.get_vector());
	EXPECT_EQ(prop2.size(),0);
}

TEST_F(AtmosphericProperty1DTest,AssignmentOperatorWorks) {
	prop6 = prop4;
	EXPECT_EQ(prop6.get_altitude_vector(),prop4.get_altitude_vector());
	EXPECT_EQ(prop6.get_vector(),prop4.get_vector());
}


TEST_F(AtmosphericProperty1DTest,AltitudeVectorRetrievable) {
	VectorWithUnits v1 = prop2.get_altitude_vector();
	EXPECT_THAT(v1, ElementsAre(
			DoubleEq(0.0),
			DoubleEq(1.0),
			DoubleEq(2.0),
			DoubleEq(3.0),
			DoubleEq(4.0),
			DoubleEq(5.0),
			DoubleEq(6.0),
			DoubleEq(7.0),
			DoubleEq(8.0),
			DoubleEq(9.0)));
}

TEST_F(AtmosphericProperty1DTest,PropertyVectorRetrievable) {
	VectorWithUnits v1 = prop2.get_vector();
	EXPECT_THAT(v1, ElementsAre(
			DoubleEq(0.0),
			DoubleEq(12.0),
			DoubleEq(8.0),
			DoubleEq(3.0),
			DoubleEq(0.0),
			DoubleEq(-1.0),
			DoubleEq(-1.3),
			DoubleEq(-2.0),
			DoubleEq(-3.0),
			DoubleEq(-5.0)));
}



TEST_F(AtmosphericProperty1DTest,UnitConversionsWork) {
	double d1 = prop2.get(2.0);
	prop2.convert_altitude_units("m");
	double d2 = prop2.get(2000.0);
	ASSERT_DOUBLE_EQ(d1,d2);
	prop2.convert_units("F");
	ASSERT_DOUBLE_EQ(prop2.get(0.0),32.0);
}

TEST_F(AtmosphericProperty1DTest,SizeFunctionIsCorrect) {
	EXPECT_EQ( prop1.size(), 0 );
	EXPECT_EQ( prop2.size(), 10 );
}

TEST_F(AtmosphericProperty1DTest,ResampleWorks) {
	prop2.resample(0.1);
	ASSERT_EQ(prop2.size(),91);
	ASSERT_DOUBLE_EQ(prop2.get(3.0),3.0);
	VectorWithUnits v1 = prop2.get_vector();
	ASSERT_DOUBLE_EQ(v1[30],3.0);
}

TEST_F(AtmosphericProperty1DTest,GetFirstDerivativeWorks) {
	EXPECT_GT(prop2.get_first_derivative(0.5),0.0);
	EXPECT_LT(prop2.get_first_derivative(1.5),0.0);
	EXPECT_NEAR(prop5.get_first_derivative(7.5),1.0,0.1);
}

TEST_F(AtmosphericProperty1DTest,GetSecondDerivativeWorks) {
	EXPECT_NEAR(prop5.get_second_derivative(5.4),0.0,0.1);
}

TEST_F(AtmosphericProperty1DTest,ResetSplinesWorks) {
	double a1 = 2.3, a2 = 7.4, a3 = 6.1;
	double d1 = prop2.get(a1), d2 = prop2.get(a2), d3 = prop2.get(a3);
	prop2.reset_splines();
	double d11 = prop2.get(a1), d21 = prop2.get(a2), d31 = prop2.get(a3);
	ASSERT_DOUBLE_EQ(d1,d11);
	ASSERT_DOUBLE_EQ(d2,d21);
	ASSERT_DOUBLE_EQ(d3,d31);
}

TEST_F(AtmosphericProperty1DTest,GetVectorAsReturnsCorrectValues) {
	double buffer[10];
	prop2.get_altitude_vector_as(buffer, "m");
	for (size_t i = 0; i < 10; i++) {
		ASSERT_DOUBLE_EQ( buffer[i], NCPA::Units::convert(alts[i],"km","m") );
	}
	prop2.get_altitude_vector_as(buffer, "F");
	for (size_t i = 0; i < 10; i++) {
		ASSERT_DOUBLE_EQ( buffer[i], NCPA::Units::convert(temps[i],"C","F") );
	}
}

/*
Example tests:
Integer Equality                EXPECT_EQ( v2.size(), 10 );
Floating Point Equality         EXPECT_DOUBLE_EQ( 40.0, 40.0 );
Vector contents                 EXPECT_THAT( v2, ElementsAre(
                                        DoubleEq(1.0),
                                        DoubleEq(2.0)) );
Vector contents all the same	EXPECT_That( v2, ElementsAre(Each(DoubleEq(1.0))) )
Exception Thrown Properly       EXPECT_THROW( {result = command();},
                                        std::out_of_range );
*/

