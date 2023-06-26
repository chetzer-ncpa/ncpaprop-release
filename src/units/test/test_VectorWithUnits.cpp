#include "NCPAUnits.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>

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
class VectorWithUnitsTest : public ::testing::Test {
 protected:
  void SetUp() override {
	  v2 = VectorWithUnits( 10, kms, KILOMETERS );
	  v3 = VectorWithUnits( 10, kms, "kilometers" );
	  ScalarWithUnits s( 10.0, CELSIUS );
	  v4 = VectorWithUnits( 10, s );
	  v5 = VectorWithUnits( 10, 10.0, CELSIUS );
	  v6 = VectorWithUnits(v4);
	  v7 = VectorWithUnits(5, svec);
  }

  // void TearDown() override {}
  VectorWithUnits v1, v2, v3, v4, v5, v6, v7;
  double kms[10] = { 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0 };
  double temps[5] = {10,20,30,40,50};
  ScalarWithUnits svec[5]= {
		  ScalarWithUnits(10.0,CELSIUS),
		  ScalarWithUnits(20.0,CELSIUS),
		  ScalarWithUnits(30.0,CELSIUS),
		  ScalarWithUnits(40.0,CELSIUS),
		  ScalarWithUnits(50.0,CELSIUS),
  };
  double all10[10] = {10,10,10,10,10,10,10,10,10,10};
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
	double buffer[10];
	v2.get_values(buffer);
	EXPECT_DOUBLE_ARRAY_EQ(10,buffer,kms);
}

TEST_F(VectorWithUnitsTest,ConstructorWithDoubleArrayAndUnitsStringWorks) {
	EXPECT_EQ( v3.size(), 10 );
	EXPECT_EQ( v3.get_units(), KILOMETERS );
	double buffer[10];
	v3.get_values(buffer);
	EXPECT_DOUBLE_ARRAY_EQ(10,buffer,kms);
}

TEST_F(VectorWithUnitsTest,ConstructorWithConstantScalarWithUnitsWorks) {
	EXPECT_THAT( v4, SizeIs(10) );
	EXPECT_EQ( v4.get_units(), CELSIUS );
	for (size_t i = 0; i < v4.size(); i++) {
		EXPECT_DOUBLE_EQ( v4[ i ].get(), 10.0 );
	}
}

TEST_F(VectorWithUnitsTest,ConstructorWithConstantDoubleWorks) {
	EXPECT_THAT( v5, SizeIs(10) );
	EXPECT_EQ( v5.get_units(), CELSIUS );
	for (size_t i = 0; i < v5.size(); i++) {
		EXPECT_DOUBLE_EQ( v5[ i ].get(), 10.0 );
	}
}

TEST_F(VectorWithUnitsTest,CopyConstructorWorks) {
	EXPECT_THAT( v6, SizeIs(10) );
	EXPECT_EQ( v6.get_units(), CELSIUS );
	for (size_t i = 0; i < v6.size(); i++) {
		EXPECT_DOUBLE_EQ( v6[ i ].get(), 10.0 );
	}
}

TEST_F(VectorWithUnitsTest,ConstructorWithScalarWithUnitsArrayWorks) {
	EXPECT_EQ( v7.size(), 5 );
	EXPECT_EQ( v7.get_units(), CELSIUS );
	double buffer[10];
	v7.get_values(buffer);
	EXPECT_DOUBLE_ARRAY_EQ(5,temps,buffer);
}

TEST_F(VectorWithUnitsTest,AssignmentOperatorWorks) {
	v3 = v7;
	EXPECT_EQ( v3.size(), v7.size() );
	EXPECT_EQ( v3.get_units(), v7.get_units() );
	for (auto v3it = v3.cbegin(), v7it = v7.cbegin();
			v3it != v3.cend() && v7it != v7.cend();
			++v3it, ++v7it) {
		EXPECT_TRUE( *v3it == *v7it );
	}
}

TEST_F(VectorWithUnitsTest,SwapWorks) {
	size_t v6size = v6.size(), v7size = v7.size();
	units_t v6u = v6.get_units(), v7u = v7.get_units();
	swap(v6,v7);

	EXPECT_EQ( v6.size(), v7size );
	EXPECT_EQ( v6.get_units(), v7u );
	double *buffer = new double[ v6.size() ];
	v6.get_values( buffer );
	EXPECT_DOUBLE_ARRAY_EQ( 5, buffer, temps );
	delete [] buffer;
	buffer = nullptr;

	EXPECT_THAT( v7, SizeIs(v6size) );
	EXPECT_EQ( v7.get_units(), v6u );
	buffer = new double[ v7.size() ];
	v7.get_values(buffer);
	EXPECT_DOUBLE_ARRAY_EQ( 10, buffer, all10 );
	delete [] buffer;
}

TEST_F(VectorWithUnitsTest, AsArrayReturnsCorrectArray) {
	ScalarWithUnits *buffer = nullptr;
	v4.as_array(buffer);
	for (auto it = v4.cbegin(); it != v4.end(); ++it) {
		EXPECT_DOUBLE_EQ( it->get(), 10.0 );
		EXPECT_EQ( it->get_units(), CELSIUS );
	}
	delete [] buffer;
}

TEST_F(VectorWithUnitsTest,ConvertUnitsCreatesCorrectValues) {
	v2.convert_units(METERS);
	for (size_t i = 0; i < 10; i++) {
		EXPECT_DOUBLE_EQ( v2[i].get(), kms[i]*1000.0 );
	}
}

TEST_F(VectorWithUnitsTest,ConvertUnitsStoresCorrectUnits) {
	v2.convert_units(METERS);
	for (size_t i = 0; i < 10; i++) {
		EXPECT_EQ(v2[i].get_units(), METERS);
	}
}

TEST_F(VectorWithUnitsTest,ConvertUnitsWithStringCreatesCorrectValues) {
	v2.convert_units("m");
	for (size_t i = 0; i < 10; i++) {
		EXPECT_DOUBLE_EQ( v2[i].get(), kms[i]*1000.0 );
	}
}

TEST_F(VectorWithUnitsTest,ConvertUnitsWithStringStoresCorrectUnits) {
	v2.convert_units("m");
	for (size_t i = 0; i < 10; i++) {
		EXPECT_EQ(v2[i].get_units(), METERS);
	}
}

TEST_F(VectorWithUnitsTest,GetUnitsReturnsCorrectUnits) {
	EXPECT_EQ( v4.get_units(), CELSIUS );
}

TEST_F(VectorWithUnitsTest,ConvertUnitsThrowsInvalidConversionCorrectly) {
	EXPECT_THROW( {v2.convert_units(CELSIUS);},
					invalid_conversion );
}

TEST_F(VectorWithUnitsTest,FillSetsValuesCorrectly) {
	size_t n = v2.size();
	v2.fill(12.0, Units::fromString("kg/m3"));
	for (size_t i = 0; i < n; i++) {
		EXPECT_DOUBLE_EQ(v2[i].get(), 12.0);
	}
}

TEST_F(VectorWithUnitsTest,FillSetsUnitsCorrectly) {
	size_t n = v2.size();
	units_t u = Units::fromString("kg/m3");
	v2.fill(12.0, u);
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(v2[i].get_units(), u);
	}
}

TEST_F(VectorWithUnitsTest,FillSetsUnitsCorrectlyFromString) {
	size_t n = v2.size();
	units_t u = Units::fromString("kg/m3");
	v2.fill(12.0, "kg/m3");
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(v2[i].get_units(), u);
	}
}

TEST_F(VectorWithUnitsTest,FillSetsValuesCorrectlyFromScalarObject) {
	size_t n = v2.size();
	units_t u = Units::fromString("kg/m3");
	ScalarWithUnits s( 12.0, u );
	v2.fill(s);
	for (size_t i = 0; i < n; i++) {
		EXPECT_DOUBLE_EQ(v2[i].get(), s.get());
	}
}

TEST_F(VectorWithUnitsTest,FillSetsUnitsCorrectlyFromScalarObject) {
	size_t n = v2.size();
	units_t u = Units::fromString("kg/m3");
	ScalarWithUnits s( 12.0, u );
	v2.fill(s);
	for (size_t i = 0; i < n; i++) {
		EXPECT_EQ(v2[i].get_units(), s.get_units());
	}
}

TEST_F(VectorWithUnitsTest,GetValuesReturnsCorrectArray) {
	double buffer[10];
	size_t n;
	v2.get_values(n,buffer);
	EXPECT_DOUBLE_ARRAY_EQ(n,buffer,kms);
}

TEST_F(VectorWithUnitsTest,GetValuesWithoutSizeReturnsCorrectArray) {
	double buffer[10];
	v2.get_values(buffer);
	EXPECT_DOUBLE_ARRAY_EQ(10,buffer,kms);
}

TEST_F(VectorWithUnitsTest,IsNormalizedReturnsTrueIfAllUnitsMatch) {
	EXPECT_TRUE( v2.is_normalized() );
}

TEST_F(VectorWithUnitsTest,IsNormalizedReturnsFalseIfNotAllUnitsMatch) {
	v2[3].set_units(CELSIUS);
	EXPECT_FALSE( v2.is_normalized() );
}

TEST_F(VectorWithUnitsTest,NormalizeUnitsCorrectlyNormalizesToUnitsOfFirstEntry) {
	v2[0].convert(METERS);
	ASSERT_FALSE( v2.is_normalized() );
	v2.normalize_units();
	ASSERT_TRUE( v2.is_normalized() );
	for (size_t i = 0; i < v2.size(); i++) {
		EXPECT_EQ(v2[i].get_units(), METERS);
	}
}

TEST_F(VectorWithUnitsTest,SetResizesCorrectly) {
	ASSERT_EQ( v2.size(), 10 );
	v2.set(5,temps,CELSIUS);
	ASSERT_EQ( v2.size(), 5 );
}

TEST_F(VectorWithUnitsTest,SetSetsNewVectorValuesCorrectly) {
	v2.set(5,temps,CELSIUS);
	for (size_t i = 0; i < 5; i++) {
		ASSERT_DOUBLE_EQ( v2[i].get(), temps[i] );
	}
}

TEST_F(VectorWithUnitsTest,SetSetsNewVectorUnitsCorrectly) {
	v2.set(5,temps,CELSIUS);
	for (size_t i = 0; i < 5; i++) {
		ASSERT_EQ( v2[i].get_units(), CELSIUS );
	}
}

TEST_F(VectorWithUnitsTest,SetSetsNewVectorUnitsCorrectlyFromString) {
	v2.set(5,temps,"C");
	for (size_t i = 0; i < 5; i++) {
		ASSERT_EQ( v2[i].get_units(), CELSIUS );
	}
}

TEST_F(VectorWithUnitsTest,SetSetsNewVectorValuesCorrectlyFromArrayOfScalarObjects) {
	v2.set(5,svec);
	for (size_t i = 0; i < 5; i++) {
		ASSERT_DOUBLE_EQ( v2[i].get(), temps[i] );
		ASSERT_EQ( v2[i].get_units(), CELSIUS );
	}
}

TEST_F(VectorWithUnitsTest,SetUnitsSetsAllUnitsCorrectly) {
	v2.set_units(METERS);
	for (size_t i = 0; i < 5; i++) {
		ASSERT_EQ( v2[i].get_units(), METERS );
	}
}

TEST_F(VectorWithUnitsTest,SetUnitsSetsAllUnitsCorrectlyFromString) {
	v2.set_units("C");
	for (size_t i = 0; i < 5; i++) {
		ASSERT_EQ( v2[i].get_units(), CELSIUS );
	}
}
