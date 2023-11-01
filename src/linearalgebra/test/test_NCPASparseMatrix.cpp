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

#ifndef EXPECT_ARRAY_EQ
#define EXPECT_ARRAY_EQ(N,A,Ex) for (int i = 0; i < N; i++) { EXPECT_EQ(A[i],Ex[i]); }
#endif

// test fixtures
class SparseMatrixTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		new_v = SparseVector<int>( {0,1,4}, {1,2,5} );

		smat1 = new SparseMatrix<int>( 5, 8 );
		smat1->set_row( 0, { 0, 4, 5, 6, 7 }, { 3, 2, 3, 2, 3 } )
			->set_row( 1, { 0, 1, 2, 3, 5, 7 }, { 3, 1, 3, 1, 2, 2 } )
			->set_row( 2, { 1, 2, 3, 4, 5, 6, 7 }, { 2, 3, 3, 3, 1, 1, 1 } )
			->set_row( 3, { 0, 1, 2, 3, 4, 5, 7 }, { 3, 3, 1, 3, 3, 2, 3 } )
			->set_row( 4, { 0, 1, 2, 3, 4 }, { 2, 3, 3, 3, 2 } );

		smat2 = new SparseMatrix<int>( 8, 5 );
		smat2->set_row( 0, { 0, 1, 2, 3, 4 }, { 1, 2, 2, 2, 2 } )
				->set_row( 1, { 0, 1, 2, 3 }, { 1, 2, 1, 2 } )
				->set_row( 2, { 0, 2, 3, 4 }, { 2, 1, 1, 2 } )
				->set_row( 3, { 0, 1, 4 }, { 2, 2, 1 } )
				->set_row( 4, { 1, 2 }, { 1, 2 } )
				->set_row( 5, { 0 }, { 1 } )
				->set_row( 6, { 0, 2, 3, 4 }, { 1, 1, 2, 1 } )
				->set_row( 7, { 0, 1, 2, 4 }, { 1, 1, 2, 1 } );

		smat2_t = new SparseMatrix<int>( 5, 8 );
		smat2_t ->set_row( 0, { 0, 1, 2, 3, 5, 6, 7 }, { 1, 1, 2, 2, 1, 1, 1 } )
				->set_row( 1, { 0, 1, 3, 4, 7 }, { 2, 2, 2, 1, 1 } )
				->set_row( 2, { 0, 1, 2, 4, 6, 7 }, { 2, 1, 1, 2, 1, 2 } )
				->set_row( 3, { 0, 1, 2, 6 }, { 2, 2, 1, 2 } )
				->set_row( 4, { 0, 2, 3, 6, 7 }, { 2, 2, 1, 1, 1 } );

		smat3 = new SparseMatrix<int>( *smat2 );
	}

	// If there's any cleanup other than normal destructors
	//void TearDown() override {}

	// class members here
	SparseMatrix<int> *smat1, *smat2, *smat2_t, *smat3, *smat4, smat5;
	Matrix<int> *mat1;
	SparseVector<int> new_v;
};

TEST_F(SparseMatrixTest,GetColumnIndicesWorksCorrectly) {
	vector<size_t> i1 = smat2->get_column_indices(4),
			expected1( { 1, 2 } );
	EXPECT_ARRAY_EQ( 2, i1, expected1 );
	vector<size_t> i2 = smat2_t->get_column_indices(3),
			expected2( { 0, 1, 2, 6 } );
	EXPECT_ARRAY_EQ( 4, i2, expected2 );
}

TEST_F(SparseMatrixTest,EqualityOperatorReturnsTrueIfEqual) {
	EXPECT_TRUE( *smat1 == *smat1 );
}

TEST_F(SparseMatrixTest,EqualityOperatorReturnsFalseIfNotEqual) {
	EXPECT_FALSE( *smat1 == *smat2 );
}

TEST_F(SparseMatrixTest,InequalityOperatorReturnsFalseIfEqual) {
	EXPECT_FALSE( *smat1 != *smat1 );
}

TEST_F(SparseMatrixTest,InequalityOperatorReturnsTrueIfNotEqual) {
	EXPECT_TRUE( *smat1 != *smat2 );
}

TEST_F(SparseMatrixTest,CopyConstructorCreatesCorrectCopy) {
	EXPECT_TRUE( *smat2 == *smat3 );
}

TEST_F(SparseMatrixTest,AssignmentOperatorCreatesCopy) {
	smat5 = *smat2;
	EXPECT_TRUE( smat5 == *smat2 );
	EXPECT_FALSE( &smat5 == smat2 );
}

TEST_F(SparseMatrixTest,CloneMethodCreatesCopy) {
	mat1 = smat2->clone();
	EXPECT_TRUE( *mat1 == *smat2 );
	EXPECT_FALSE( mat1 == smat2 );
	delete mat1;
}

// Tests that use the fixture
TEST_F(SparseMatrixTest,MatricesHaveCorrectRowCounts) {
	ASSERT_EQ( smat1->rows(), 5 );
	ASSERT_EQ( smat2->rows(), 8 );
}

TEST_F(SparseMatrixTest,MatricesHaveCorrectColumnCounts) {
	ASSERT_EQ( smat1->columns(), 8 );
	ASSERT_EQ( smat2->columns(), 5 );
}

TEST_F(SparseMatrixTest,InitializeSetsRowsCorrectly) {
	smat1->initialize( 3, 3 );
	EXPECT_EQ( smat1->rows(), 3 );
}

TEST_F(SparseMatrixTest,InitializeSetsColumnsCorrectly) {
	smat1->initialize( 3, 3 );
	EXPECT_EQ( smat1->columns(), 3 );
}

TEST_F(SparseMatrixTest,GetRetrievesCorrectValues) {
	EXPECT_EQ( smat1->get( 2, 5 ), 1 );
	EXPECT_EQ( smat2->get( 0, 3 ), 2 );
}

TEST_F(SparseMatrixTest,SetAltersValuesCorrectly) {
	ASSERT_EQ( smat1->get( 2, 5 ), 1 );
	smat1->set( 2, 5, 10 );
	ASSERT_EQ( smat1->get( 2, 5 ), 10 );
}

TEST_F(SparseMatrixTest,GetRowGetsRowCorrectlyAsArray) {
	int *row = zeros<int>( 8 );
	size_t n = 8;
	int expected[8] = { 0, 2, 3, 3, 3, 1, 1, 1 };
	smat1->get_row( 2, n, row );
	EXPECT_ARRAY_EQ( 8, row, expected );
	delete [] row;
}

TEST_F(SparseMatrixTest,GetRowCreatesArrayIfPassedNullptr) {
	int *row = nullptr;
	size_t n = 0;
	int expected[8] = { 0, 2, 3, 3, 3, 1, 1, 1 };
	smat1->get_row( 2, n, row );
	EXPECT_EQ( n, 8 );
	EXPECT_ARRAY_EQ( 8, row, expected );
	delete [] row;
}

TEST_F(SparseMatrixTest,GetRowThrowsOutOfRangeIfWrongSizeArrayProvided) {
	int *row = zeros<int>( 6 );
	size_t n = 6;
	EXPECT_THROW( {smat1->get_row( 2, n, row );}, out_of_range );
}

TEST_F(SparseMatrixTest,GetRowGetsRowCorrectlyAsVector) {
	vector<int> row;
	int expected[8] = { 0, 2, 3, 3, 3, 1, 1, 1 };
	smat1->get_row( 2, row );
	EXPECT_EQ( row.size(), 8 );
	EXPECT_ARRAY_EQ( 8, row, expected );
}

TEST_F(SparseMatrixTest,GetRowReturnsCorrectSparseVector) {
	SparseVector<int> v = smat1->get_row( 2 );
	EXPECT_EQ( v.size(), 7 );
	int expected[8] = { 0, 2, 3, 3, 3, 1, 1, 1 };
	for (auto vit = v.cbegin(); vit != v.cend(); ++vit) {
		EXPECT_EQ( vit->second, expected[vit->first] );
	}
}

TEST_F(SparseMatrixTest,GetRowReturnsCorrectSparseVectorAsArgument) {
	SparseVector<int> v;
	smat1->get_row( 2, v );
	EXPECT_EQ( v.size(), 7 );
	int expected[8] = { 0, 2, 3, 3, 3, 1, 1, 1 };
	for (auto vit = v.cbegin(); vit != v.cend(); ++vit) {
		EXPECT_EQ( vit->second, expected[vit->first] );
	}
}


TEST_F(SparseMatrixTest,GetColumnGetsColumnCorrectlyAsArray) {
	int *col = zeros<int>( 5 );
	size_t n = 5;
	int expected[5] = { 3, 3, 0, 3, 2 };
	smat1->get_column( 0, n, col );
	EXPECT_ARRAY_EQ( 5, col, expected );
	delete [] col;
}

TEST_F(SparseMatrixTest,GetColumnCreatesArrayIfPassedNullptr) {
	int *col = nullptr;
	size_t n = 0;
	int expected[5] = { 3, 3, 0, 3, 2 };
	smat1->get_column( 0, n, col );
	EXPECT_EQ( n, 5 );
	EXPECT_ARRAY_EQ( 5, col, expected );
	delete [] col;
}

TEST_F(SparseMatrixTest,GetColumnThrowsOutOfRangeIfWrongSizeArrayProvided) {
	int *col = zeros<int>( 6 );
	size_t n = 6;
	EXPECT_THROW( {smat1->get_column( 2, n, col );}, out_of_range );
}

TEST_F(SparseMatrixTest,GetColumnGetsColumnCorrectlyAsVector) {
	vector<int> col;
	int expected[5] = { 3, 3, 0, 3, 2 };
	smat1->get_column( 0, col );
	EXPECT_EQ( col.size(), 5 );
	EXPECT_ARRAY_EQ( 5, col, expected );
}

TEST_F(SparseMatrixTest,GetColumnReturnsCorrectSparseVector) {
	SparseVector<int> v = smat1->get_column( 0 );
	EXPECT_EQ( v.size(), 4 );
	int expected[5] = { 3, 3, 0, 3, 2 };
	for (auto vit = v.cbegin(); vit != v.cend(); ++vit) {
		EXPECT_EQ( vit->second, expected[vit->first] );
	}
}

TEST_F(SparseMatrixTest,GetColumnReturnsCorrectSparseVectorAsArgument) {
	SparseVector<int> v;
	smat1->get_column( 0, v );
	EXPECT_EQ( v.size(), 4 );
	int expected[5] = { 3, 3, 0, 3, 2 };
	for (auto vit = v.cbegin(); vit != v.cend(); ++vit) {
		EXPECT_EQ( vit->second, expected[vit->first] );
	}
}



TEST_F(SparseMatrixTest,SetRowWithArraySetsRowCorrectly) {
	size_t *inds = NCPA::zeros<size_t>( 5 );
	int *vals = NCPA::zeros<int>( 5 );
	for (size_t i = 0; i < 5; i++) {
		inds[i] = i;
		vals[i] = i+1;
	}
	smat2->set_row( 0, 5, inds, vals );

	vector<int> row0;
	smat2->get_row( 0, row0 );
	EXPECT_ARRAY_EQ( 5, row0, vals );
}

TEST_F(SparseMatrixTest,SetRowReturnsThis) {
	size_t *inds = NCPA::zeros<size_t>( 5 );
	int *vals = NCPA::zeros<int>( 5 );
	for (size_t i = 0; i < 5; i++) {
		inds[i] = i;
		vals[i] = i+1;
	}
	NCPA::SparseMatrix<int> *ptr = static_cast<NCPA::SparseMatrix<int>*>(
			smat2->set_row( 0, 5, inds, vals ) );
	EXPECT_EQ(ptr,smat2);
}

TEST_F(SparseMatrixTest,SetRowWithVectorSetsRowCorrectly) {
	vector<int> vals(5);
	for (size_t i = 0; i < 5; i++) {
		vals[i] = i+1;
	}
	smat2->set_row( 0, vals );

	vector<int> row0;
	smat2->get_row( 0, row0 );
	EXPECT_ARRAY_EQ( 5, row0, vals );
}

TEST_F(SparseMatrixTest,SetRowWithTwoVectorsSetsRowCorrectly) {
	vector<int> vals = { 1, 2, 5 };
	vector<size_t> inds = {0, 1, 3 };
	smat2->set_row( 0, inds, vals );

	vector<int> row0;
	smat2->get_row( 0, row0 );
	vector<int> expected = { 1, 2, 0, 5, 0 };
	EXPECT_ARRAY_EQ( 5, row0, expected );
}

TEST_F(SparseMatrixTest,SetRowWithSparseVectorSetsRowCorrectly) {
	smat2->set_row( 0, new_v );
	vector<int> row0;
	smat2->get_row( 0, row0 );
	vector<int> expected = { 1, 2, 0, 0, 5 };
	EXPECT_ARRAY_EQ( 5, row0, expected );
}

TEST_F(SparseMatrixTest,SetColumnWithArraySetsColumnCorrectly) {
	size_t *inds = NCPA::zeros<size_t>( 5 );
	int *vals = NCPA::zeros<int>( 5 );
	for (size_t i = 0; i < 5; i++) {
		inds[i] = i;
		vals[i] = i+1;
	}
	smat1->set_column( 0, 5, inds, vals );

	vector<int> col0;
	smat1->get_column( 0, col0 );
	EXPECT_ARRAY_EQ( 5, col0, vals );
}

TEST_F(SparseMatrixTest,SetColumnReturnsThis) {
	size_t *inds = NCPA::zeros<size_t>( 5 );
	int *vals = NCPA::zeros<int>( 5 );
	for (size_t i = 0; i < 5; i++) {
		inds[i] = i;
		vals[i] = i+1;
	}
	NCPA::SparseMatrix<int> *ptr = static_cast<NCPA::SparseMatrix<int>*>(
			smat1->set_column( 0, 5, inds, vals ) );
	EXPECT_EQ(ptr,smat1);
}

TEST_F(SparseMatrixTest,SetColumnWithVectorSetsColumnCorrectly) {
	vector<int> vals(5);
	for (size_t i = 0; i < 5; i++) {
		vals[i] = i+1;
	}
	smat1->set_column( 0, vals );

	vector<int> col0;
	smat1->get_column( 0, col0 );
	EXPECT_ARRAY_EQ( 5, col0, vals );
}

TEST_F(SparseMatrixTest,SetColumnWithTwoVectorsSetsColumnCorrectly) {
	vector<int> vals = { 1, 2, 5 };
	vector<size_t> inds = {0, 1, 3 };
	smat1->set_column( 0, inds, vals );

	vector<int> col0;
	smat1->get_column( 0, col0 );
	vector<int> expected = { 1, 2, 0, 5, 0 };
	EXPECT_ARRAY_EQ( 5, col0, expected );
}

TEST_F(SparseMatrixTest,SetColumnWithSparseVectorSetsColumnCorrectly) {
	smat1->set_column( 0, new_v );
	vector<int> col0;
	smat1->get_column( 0, col0 );
	vector<int> expected = { 1, 2, 0, 0, 5 };
	EXPECT_ARRAY_EQ( 5, col0, expected );
}

TEST_F(SparseMatrixTest,IsReadyReturnsTrueForPopulatedMatrix) {
	EXPECT_TRUE( smat1->is_ready() );
}

TEST_F(SparseMatrixTest,IsReadyReturnsFalseForUnpopulatedMatrix) {
	EXPECT_FALSE( smat5.is_ready() );
}

TEST_F(SparseMatrixTest,AdditionWorksCorrectly) {
	smat2->add( smat3 );
	for (size_t row = 0; row < smat2->rows(); row++) {
		SparseVector<int> v = smat2->get_row(row);
		for (auto it = v.begin(); it != v.end(); ++it) {
			EXPECT_EQ( it->second, 2*smat3->get(row,it->first) );
		}
	}
}

TEST_F(SparseMatrixTest,AdditionPreservesStorageEfficiency) {
	smat2->add( smat3 );
	for (size_t row = 0; row < smat2->rows(); row++) {
		SparseVector<int> v2 = smat2->get_row(row),
				v3 = smat3->get_row(row);
		EXPECT_EQ( v2.size(), v3.size() );
	}
}

TEST_F(SparseMatrixTest,AdditionThrowsOutOfRangeIfSizesDoNotMatch) {
	EXPECT_THROW( { smat2->add( smat1 ); }, out_of_range );
}

TEST_F(SparseMatrixTest,SubtractionWorksCorrectly) {
	smat2->subtract( smat3 );
	for (size_t row = 0; row < smat2->rows(); row++) {
		SparseVector<int> v = smat2->get_row(row);
		for (auto it = v.begin(); it != v.end(); ++it) {
			EXPECT_EQ( it->second, 0 );
		}
	}
}

TEST_F(SparseMatrixTest,SubtractionThrowsOutOfRangeIfSizesDoNotMatch) {
	EXPECT_THROW( { smat2->subtract( smat1 ); }, out_of_range );
}

TEST_F(SparseMatrixTest,TransposeWorksCorrectly) {
//	cout << "Before transpose: ";
//	smat2->print();

	smat2->transpose();
//	cout << "After transpose: ";
//	smat2->print();
//	cout << "smat2_t: ";
//	smat2_t->print();
	EXPECT_TRUE( *smat2 == *smat2_t );
}

TEST_F(SparseMatrixTest,ScaleWorksCorrectly) {
	smat3->add(smat2)->add(smat2);
	smat2->scale(3.0);
	EXPECT_TRUE( *smat2 == *smat3 );
}

TEST_F(SparseMatrixTest,ScaleWorksCorrectlyForNegativeScalars) {
	smat3->subtract(smat2)->subtract(smat2);
	smat2->scale(-1.0);
	EXPECT_TRUE( *smat2 == *smat3 );
}

// Tests that don't use the fixture
//TEST(SparseMatrixTest,TestName2) {
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
