#include "NCPALinearAlgebra.h"
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>

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
class SparseVectorTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		v1 = SparseVector<int>(5);
		v1[0] = 1;
		v1[2] = 3;
		v1[99] = 100;
		v1[3] = 4;
		v1[4] = 5;
	}

	// If there's any cleanup other than normal destructors
	//void TearDown() override {}

	// class members here
	SparseVector<int> v1, v2;
};


// Tests that use the fixture
TEST_F(SparseVectorTest,VectorReportsCorrectSize) {
	EXPECT_EQ( v1.size(), 5 );
	EXPECT_EQ( v2.size(), 0 );
}

TEST_F(SparseVectorTest,VectorRetrievesCorrectValues) {
	EXPECT_EQ( v1[2], 3 );
	EXPECT_EQ( v1[99], 100 );
	EXPECT_EQ( v1.at( 2 ), 3 );
	EXPECT_EQ( v1.at( 99 ), 100 );
}

TEST_F(SparseVectorTest,VectorAddsZeroValueWhenCalledWithBrackets) {
	EXPECT_EQ( v1.size(), 5 );
	EXPECT_EQ( v1[1], 0 );
	EXPECT_EQ( v1.size(), 6 );
}

TEST_F(SparseVectorTest,VectorThrowsOutOfRangeWhenAtCalledForMissingValue) {
	EXPECT_THROW( {int i = v1.at(1); }, out_of_range );
}

TEST_F(SparseVectorTest,ConstructorWorksWithVectors) {
	vector<size_t> inds( { 0, 2, 3, 4, 99 } );
	vector<int> vals( { 1, 3, 4, 5, 100 } );
	SparseVector<int> v( inds, vals );
	EXPECT_EQ( v1, v );
}

TEST_F(SparseVectorTest,ConstructorWorksWithArrays) {
	size_t inds[] = { 0, 2, 3, 4, 99 };
	int vals[] = { 1, 3, 4, 5, 100 };
	SparseVector<int> v( 5, inds, vals );
	EXPECT_EQ( v1, v );
}

TEST_F(SparseVectorTest,ConstructorWorksWithInitializerLists) {
	SparseVector<int> v(
			{ 0, 2, 3, 4, 99 },
			{ 1, 3, 4, 5, 100 } );
	EXPECT_EQ( v1, v );
}

TEST_F(SparseVectorTest,ConstructorWorksWithVectorAndConstant) {
	SparseVector<int> v( vector<size_t>( { 0, 2, 3, 4, 99 } ), 1 );
	EXPECT_EQ( v[0], 1 );
	EXPECT_EQ( v[2], 1 );
	EXPECT_EQ( v[3], 1 );
	EXPECT_EQ( v[4], 1 );
	EXPECT_EQ( v[99], 1 );

}

TEST_F(SparseVectorTest,ConstructorWorksWithArrayAndConstant) {
	size_t inds[] = { 0, 2, 3, 4, 99 };
	SparseVector<int> v( 5, inds, 1 );
	EXPECT_EQ( v[0], 1 );
	EXPECT_EQ( v[2], 1 );
	EXPECT_EQ( v[3], 1 );
	EXPECT_EQ( v[4], 1 );
	EXPECT_EQ( v[99], 1 );
}

TEST_F(SparseVectorTest,ConstructorWorksWithInitializerListAndConstant) {
	SparseVector<int> v( { 0, 2, 3, 4, 99 }, 1 );
	EXPECT_EQ( v[0], 1 );
	EXPECT_EQ( v[2], 1 );
	EXPECT_EQ( v[3], 1 );
	EXPECT_EQ( v[4], 1 );
	EXPECT_EQ( v[99], 1 );
}

// Tests that don't use the fixture
//TEST(SparseVectorTest,TestName2) {
//
//}

/*
Example tests:
Integer Equality		EXPECT_EQ( v2.size(), 10 );
Floating Point Equality		EXPECT_DOUBLE_EQ( 40.0, 40.0 );
Floating Point Rounding		EXPECT_NEAR( 40.123, 40.12, 1e-2 );
Vector contents			EXPECT_THAT( v2, ElementsAre(
					DoubleEq(1.0),
					DoubleEq(2.0)) );
Exception Thrown Properly	EXPECT_THROW( {result = command();},
					std::out_of_range );
*/
