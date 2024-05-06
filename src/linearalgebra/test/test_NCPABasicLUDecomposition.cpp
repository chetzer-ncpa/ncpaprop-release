#include "NCPABasicLUDecomposition.h"
#include "NCPAMatrix.h"
#include <gmock/gmock.h>
#include <cmath>
#include <cstdlib>
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
class BasicLUDecompositionTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		dmat = new BasicMatrix<double>( 4, 4 );
		for (size_t r = 0; r < 4; r++) {
			for (size_t c = 0; c < 4; c++) {
				dmat->set( r, c, (double)(rand() % 10) );
			}
		}
		dmat->finalize();

		decomp = new BasicLUDecomposition<double>();
	}

	// If there's any cleanup other than normal destructors
	//void TearDown() override {}

	// class members here
	BasicLUDecomposition<double> *decomp;
	BasicMatrix<double> *dmat, *lower, *upper, *permut;
};

TEST_F(BasicLUDecompositionTest,DecompositionIsNotReadyOnBuild) {
	EXPECT_FALSE( decomp->is_ready() );
}

TEST_F(BasicLUDecompositionTest,InitializeMakesDecompositionReady) {
	decomp->initialize( 4, 4 );
	EXPECT_TRUE( decomp->is_ready() );
}

TEST_F(BasicLUDecompositionTest,InitializeCreatesLowerMatrix) {
	decomp->initialize( 4, 4 );
	EXPECT_EQ( decomp->lower()->rows(), 4 );
	EXPECT_EQ( decomp->lower()->columns(), 4 );
}

TEST_F(BasicLUDecompositionTest,InitializeCreatesUpperMatrix) {
	decomp->initialize( 4, 4 );
	EXPECT_EQ( decomp->upper()->rows(), 4 );
	EXPECT_EQ( decomp->upper()->columns(), 4 );
}

TEST_F(BasicLUDecompositionTest,InitializeCreatesPermutationMatrix) {
	decomp->initialize( 4, 4 );
	EXPECT_EQ( decomp->permutation()->rows(), 4 );
	EXPECT_EQ( decomp->permutation()->columns(), 4 );
}

TEST_F(BasicLUDecompositionTest,CalculateWorksCorrectly) {
	decomp->initialize( 4, 4 );
	decomp->compute( dmat );
	NCPA::Matrix<double> *left = new NCPA::BasicMatrix<double>(),
						 *right = new NCPA::BasicMatrix<double>();
	decomp->permutation()->multiply( dmat, left );
	decomp->lower()->multiply( decomp->upper(), right );
	for (size_t r = 0; r < left->rows(); r++) {
		for (size_t c = 0; c < left->columns(); c++) {
			EXPECT_NEAR( left->get(r,c), right->get(r,c), 1e-10 );
		}
	}
//	cout << "Lower with pivot:" << endl;
//	decomp->lower()->print();
//	cout << "Upper with pivot:" << endl;
//	decomp->upper()->print();
//	cout << "Permutation with pivot:" << endl;
//	decomp->permutation()->print();
}

TEST_F(BasicLUDecompositionTest,CalculateWorksCorrectlyWithoutPivot) {
	decomp->initialize( 4, 4 );
	decomp->compute( dmat, false );
	NCPA::Matrix<double> *left = new NCPA::BasicMatrix<double>(),
						 *right = new NCPA::BasicMatrix<double>();
	decomp->permutation()->multiply( dmat, left );
	decomp->lower()->multiply( decomp->upper(), right );
	for (size_t r = 0; r < left->rows(); r++) {
		for (size_t c = 0; c < left->columns(); c++) {
			EXPECT_NEAR( left->get(r,c), right->get(r,c), 1e-10 );
		}
	}

	// permutation should be the identity matrix
	EXPECT_TRUE( decomp->permutation()->is_diagonal() );
//	cout << "Lower without pivot:" << endl;
//	decomp->lower()->print();
//	cout << "Upper without pivot:" << endl;
//	decomp->upper()->print();
//	cout << "Permutation without pivot:" << endl;
//	decomp->permutation()->print();
}

TEST_F(BasicLUDecompositionTest,ClearSetsLowerToNull) {
	decomp->initialize(4,4);
	EXPECT_NE( decomp->lower(), nullptr );
	decomp->clear();
	EXPECT_EQ( decomp->lower(), nullptr );
}

TEST_F(BasicLUDecompositionTest,ClearSetsUpperToNull) {
	decomp->initialize(4,4);
	EXPECT_NE( decomp->upper(), nullptr );
	decomp->clear();
	EXPECT_EQ( decomp->upper(), nullptr );
}

TEST_F(BasicLUDecompositionTest,ClearSetsPermutationToNull) {
	decomp->initialize(4,4);
	EXPECT_NE( decomp->permutation(), nullptr );
	decomp->clear();
	EXPECT_EQ( decomp->permutation(), nullptr );
}

// Tests that use the fixture

// Tests that don't use the fixture
//TEST(BasicLUDecompositionTest,TestName2) {
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
