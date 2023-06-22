#include "NCPAUnits.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace NCPA;
using namespace testing;


// test fixture
class VectorWithUnitsTest : public ::testing::Test {
 protected:
  void SetUp() override {
	  v2 = VectorWithUnits( 10, kms, KILOMETERS );
	  v3 = VectorWithUnits( 10, kms, "kilometers" );
	  ScalarWithUnits s( 10.0, CELSIUS );
	  v4 = VectorWithUnits( 10, s );
	  v5 = VectorWithUnits( 10, 10.0, CELSIUS );
	  v6 = VectorWithUnits(v4);
	  ScalarWithUnits svec[5] = {
			  ScalarWithUnits(10.0,CELSIUS),
			  ScalarWithUnits(20.0,CELSIUS),
			  ScalarWithUnits(30.0,CELSIUS),
			  ScalarWithUnits(40.0,CELSIUS),
			  ScalarWithUnits(50.0,CELSIUS),
	  };
	  v7 = VectorWithUnits(5, svec);
  }

  // void TearDown() override {}
  VectorWithUnits v1, v2, v3, v4, v5, v6, v7;
  double kms[10] = { 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0 };
  const units_t KILOMETERS 	= Units::fromString("km"),
		  	    METERS 		= Units::fromString("m"),
				CELSIUS 	= Units::fromString("C");
};




TEST_F(VectorWithUnitsTest,DefaultConstructorIsEmpty) {
	ASSERT_EQ( v1.size(), 0 );
}

TEST_F(VectorWithUnitsTest,ConstructorWithDoubleArrayAndUnitsObjectWorks) {
	EXPECT_EQ( v2.size(), 10 );
	EXPECT_EQ( v2.get_units(), KILOMETERS );
	EXPECT_THAT( v2, ElementsAre(
			DoubleEq(1.0),
			DoubleEq(2.0),
			DoubleEq(3.0),
			DoubleEq(4.0),
			DoubleEq(5.0),
			DoubleEq(4.0),
			DoubleEq(3.0),
			DoubleEq(2.0),
			DoubleEq(1.0),
			DoubleEq(0.0)) );
}

TEST_F(VectorWithUnitsTest,ConstructorWithDoubleArrayAndUnitsStringWorks) {
	EXPECT_EQ( v3.size(), 10 );
	EXPECT_EQ( v3.get_units(), KILOMETERS );
	EXPECT_THAT( v3, ElementsAre(
				DoubleEq(1.0),
				DoubleEq(2.0),
				DoubleEq(3.0),
				DoubleEq(4.0),
				DoubleEq(5.0),
				DoubleEq(4.0),
				DoubleEq(3.0),
				DoubleEq(2.0),
				DoubleEq(1.0),
				DoubleEq(0.0)) );
}

TEST_F(VectorWithUnitsTest,ConstructorWithConstantScalarWithUnitsWorks) {
	EXPECT_THAT( v4, SizeIs(10) );
	EXPECT_EQ( v4.get_units(), CELSIUS );
	EXPECT_THAT( v4, Each(DoubleEq(10.0)) );
}

TEST_F(VectorWithUnitsTest,ConstructorWithConstantDoubleWorks) {
	EXPECT_THAT( v5, SizeIs(10) );
	EXPECT_EQ( v5.get_units(), CELSIUS );
	EXPECT_THAT( v5, Each(DoubleEq(10.0)) );
}

TEST_F(VectorWithUnitsTest,CopyConstructorWorks) {
	EXPECT_THAT( v6, SizeIs(10) );
	EXPECT_EQ( v6.get_units(), CELSIUS );
	EXPECT_THAT( v6, Each(DoubleEq(10.0)) );
}

TEST_F(VectorWithUnitsTest,ConstructorWithScalarWithUnitsArrayWorks) {
	EXPECT_EQ( v7.size(), 5 );
	EXPECT_EQ( v7.get_units(), CELSIUS );
	EXPECT_THAT( v7, ElementsAre(
				DoubleEq(10.0),
				DoubleEq(20.0),
				DoubleEq(30.0),
				DoubleEq(40.0),
				DoubleEq(50.0)) );
}

TEST_F(VectorWithUnitsTest,AssignmentOperatorWorks) {
	v3 = v7;
	EXPECT_EQ( v3.size(), 5 );
	EXPECT_EQ( v3.get_units(), CELSIUS );
	EXPECT_THAT( v3, ElementsAre(
				DoubleEq(10.0),
				DoubleEq(20.0),
				DoubleEq(30.0),
				DoubleEq(40.0),
				DoubleEq(50.0)) );
}

TEST_F(VectorWithUnitsTest,SwapWorks) {
	swap(v6,v7);
	EXPECT_EQ( v6.size(), 5 );
	EXPECT_EQ( v6.get_units(), CELSIUS );
	EXPECT_THAT( v6, ElementsAre(
				DoubleEq(10.0),
				DoubleEq(20.0),
				DoubleEq(30.0),
				DoubleEq(40.0),
				DoubleEq(50.0)) );
	EXPECT_THAT( v7, SizeIs(10) );
	EXPECT_EQ( v7.get_units(), CELSIUS );
	EXPECT_THAT( v7, Each(DoubleEq(10.0)) );
}

TEST_F(VectorWithUnitsTest,GetUnitsWorks) {
	EXPECT_EQ( v4.get_units(), CELSIUS );
}

TEST_F(VectorWithUnitsTest,ConvertUnitsWorks) {
	v2.convert_units(METERS);
	EXPECT_THAT( v2, ElementsAre(
					DoubleEq(1000.0),
					DoubleEq(2000.0),
					DoubleEq(3000.0),
					DoubleEq(4000.0),
					DoubleEq(5000.0),
					DoubleEq(4000.0),
					DoubleEq(3000.0),
					DoubleEq(2000.0),
					DoubleEq(1000.0),
					DoubleEq(0.0)) );
	v5.convert_units("F");
	EXPECT_THAT( v5, Each(DoubleEq(50.0)) );
}

TEST_F(VectorWithUnitsTest,ConvertUnitsThrowsInvalidConversion) {
	EXPECT_THROW( {v2.convert_units(CELSIUS);},
					invalid_conversion );
}

TEST_F(VectorWithUnitsTest, SetWorks) {
	v1.set(10,kms,KILOMETERS);
	EXPECT_EQ( v1.size(), 10 );
	EXPECT_EQ( v1.get_units(), KILOMETERS );
	EXPECT_THAT( v1, ElementsAre(
			DoubleEq(1.0),
			DoubleEq(2.0),
			DoubleEq(3.0),
			DoubleEq(4.0),
			DoubleEq(5.0),
			DoubleEq(4.0),
			DoubleEq(3.0),
			DoubleEq(2.0),
			DoubleEq(1.0),
			DoubleEq(0.0)) );
}

TEST_F(VectorWithUnitsTest,SetUnitsWorks) {
	v2.set_units(METERS);
	EXPECT_EQ( v2.get_units(), METERS );
	EXPECT_THAT( v2, ElementsAre(
				DoubleEq(1.0),
				DoubleEq(2.0),
				DoubleEq(3.0),
				DoubleEq(4.0),
				DoubleEq(5.0),
				DoubleEq(4.0),
				DoubleEq(3.0),
				DoubleEq(2.0),
				DoubleEq(1.0),
				DoubleEq(0.0)) );
}

TEST_F(VectorWithUnitsTest,FillWorks) {
	v5.fill(26.0);
	EXPECT_THAT( v5, SizeIs(10) );
	EXPECT_EQ( v5.get_units(), CELSIUS );
	EXPECT_THAT( v5, Each(DoubleEq(26.0)) );
}

TEST_F(VectorWithUnitsTest,FillWithNewSizeWorks) {
	v6.fill(5,12.0);
	EXPECT_THAT( v6, SizeIs(5) );
	EXPECT_EQ( v6.get_units(), CELSIUS );
	EXPECT_THAT( v6, Each(DoubleEq(12.0)) );
}


TEST_F(VectorWithUnitsTest,SetValuesWorks) {
	v7.set_values( 10, kms );
	EXPECT_EQ( v7.size(), 10 );
	EXPECT_EQ( v7.get_units(), CELSIUS );
	EXPECT_THAT( v7, ElementsAre(
				DoubleEq(1.0),
				DoubleEq(2.0),
				DoubleEq(3.0),
				DoubleEq(4.0),
				DoubleEq(5.0),
				DoubleEq(4.0),
				DoubleEq(3.0),
				DoubleEq(2.0),
				DoubleEq(1.0),
				DoubleEq(0.0)) );
}

TEST_F(VectorWithUnitsTest,GetVectorWorks) {
	double *km_copy = new double[10];
	units_t km_units;
	v2.get_vector(km_copy, km_units);
	EXPECT_DOUBLE_EQ(km_copy[0],1.0);
	EXPECT_DOUBLE_EQ(km_copy[1],2.0);
	EXPECT_DOUBLE_EQ(km_copy[2],3.0);
	EXPECT_DOUBLE_EQ(km_copy[3],4.0);
	EXPECT_DOUBLE_EQ(km_copy[4],5.0);
	EXPECT_DOUBLE_EQ(km_copy[5],4.0);
	EXPECT_DOUBLE_EQ(km_copy[6],3.0);
	EXPECT_DOUBLE_EQ(km_copy[7],2.0);
	EXPECT_DOUBLE_EQ(km_copy[8],1.0);
	EXPECT_DOUBLE_EQ(km_copy[9],0.0);
	EXPECT_EQ( km_units, KILOMETERS );

	v7.get_vector(km_copy);
	EXPECT_DOUBLE_EQ(km_copy[0],10.0);
	EXPECT_DOUBLE_EQ(km_copy[1],20.0);
	EXPECT_DOUBLE_EQ(km_copy[2],30.0);
	EXPECT_DOUBLE_EQ(km_copy[3],40.0);
	EXPECT_DOUBLE_EQ(km_copy[4],50.0);
	EXPECT_DOUBLE_EQ(km_copy[5],4.0);
	EXPECT_DOUBLE_EQ(km_copy[6],3.0);
	EXPECT_DOUBLE_EQ(km_copy[7],2.0);
	EXPECT_DOUBLE_EQ(km_copy[8],1.0);
	EXPECT_DOUBLE_EQ(km_copy[9],0.0);
	delete [] km_copy;
}
