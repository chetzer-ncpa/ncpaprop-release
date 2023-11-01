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
class MatrixTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
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

		smat3 = new SparseMatrix<int>( *smat2 );

		nsmat1 = new SparseMatrix<int>( 5, 8 );
		nsmat1->set_row( 0, { 0, 4, 5, 6, 7 }, { -3, -2, -3, -2, -3 } )
			->set_row( 1, { 0, 1, 2, 3, 5, 7 }, { -3, -1, -3, -1, -2, -2 } )
			->set_row( 2, { 1, 2, 3, 4, 5, 6, 7 }, { -2, -3, -3, -3, -1, -1, -1 } )
			->set_row( 3, { 0, 1, 2, 3, 4, 5, 7 }, { -3, -3, -1, -3, -3, -2, -3 } )
			->set_row( 4, { 0, 1, 2, 3, 4 }, { -2, -3, -3, -3, -2 } );

		mat1 = static_cast<Matrix<int>*>( smat1 );
		mat2 = static_cast<Matrix<int>*>( smat2 );
		mat3 = static_cast<Matrix<int>*>( smat3 );
		nmat1 = static_cast<Matrix<int>*>( nsmat1 );
	}

	// If there's any cleanup other than normal destructors
	//void TearDown() override {}

	// class members here
	SparseMatrix<int> *smat1, *smat2, *smat3, *smat4, *nsmat1;
	Matrix<int> *mat1, *mat2, *mat3, *mat4 = nullptr, *nmat1;
};

// Tests that use the fixture
TEST_F(MatrixTest,MatrixEqualityReturnsTrueIfEqual) {
	EXPECT_TRUE( *smat3 == *smat2 );
}

TEST_F(MatrixTest,MatrixEqualityReturnsFalseIfNotEqual) {
	EXPECT_FALSE( *smat2 == *smat1 );
}

TEST_F(MatrixTest,MatrixInequalityReturnsFalseIfEqual) {
	EXPECT_FALSE( *smat3 != *smat2 );
}

TEST_F(MatrixTest,MatrixEqualityReturnsTrueIfNotEqual) {
	EXPECT_TRUE( *smat2 != *smat1 );
}

TEST_F(MatrixTest,CheckDimensionsThrowsNothingForDimensionsInRange) {
	EXPECT_NO_THROW( { smat1->check_dimensions(1,1); } );
}

TEST_F(MatrixTest,CheckDimensionsThrowsOutOfRangeForRowNotInRange) {
	EXPECT_THROW( { smat1->check_dimensions(12,1); }, out_of_range );
}

TEST_F(MatrixTest,CheckDimensionsThrowsOutOfRangeForColumnNotInRange) {
	EXPECT_THROW( { smat1->check_dimensions(1,12); }, out_of_range );
}

TEST_F(MatrixTest,MultiplyThrowsOutOfRangeIfInsideDimensionsDoNotMatch) {
	EXPECT_THROW( { Matrix<int>::multiply( mat1, mat1, mat4 );}, out_of_range );
}

TEST_F(MatrixTest,MultiplyThrowsLogicErrorIfProductMatrixIsNull) {
	EXPECT_THROW( { Matrix<int>::multiply( mat1, mat2, mat4 );}, logic_error );
}

TEST_F(MatrixTest,MultiplyGivesCorrectProduct) {
	mat4 = new SparseMatrix<int>();
	Matrix<int>::multiply( mat1, mat2, mat4 );
	std::vector<int> testrow;
	double expected0[5] = { 11, 11, 18, 10, 11 },
			expected1[5] = { 16, 12, 14, 11, 15 },
			expected2[5] = { 17, 14, 14, 9, 11 },
			expected3[5] = { 19, 24, 22, 13, 14 },
			expected4[5] = { 17, 18, 14, 13, 13 };
	mat4->get_row( 0, testrow );
	EXPECT_ARRAY_EQ( 5, testrow, expected0 );
	mat4->get_row( 1, testrow );
	EXPECT_ARRAY_EQ( 5, testrow, expected1 );
	mat4->get_row( 2, testrow );
	EXPECT_ARRAY_EQ( 5, testrow, expected2 );
	mat4->get_row( 3, testrow );
	EXPECT_ARRAY_EQ( 5, testrow, expected3 );
	mat4->get_row( 4, testrow );
	EXPECT_ARRAY_EQ( 5, testrow, expected4 );
	delete mat4;
}

TEST_F(MatrixTest,NegativeWorksProperly) {
	mat1->negative();
	EXPECT_TRUE( *mat1 == *nmat1 );
}

TEST_F(MatrixTest,GetDiagonalGetsCorrectVector) {
	vector<int> diag = mat1->get_diagonal();
	int expected1[5] = { 3, 1, 3, 3, 2 };
	EXPECT_ARRAY_EQ( 5, diag, expected1 );
	diag = mat2->get_diagonal();
	int expected2[5] = { 1, 2, 1, 0, 0 };
	EXPECT_ARRAY_EQ( 5, diag, expected2 );
}

TEST_F(MatrixTest,GetDiagonalGetsCorrectArray) {
	int *diag = NCPA::zeros<int>(5);
	size_t n = 0;
	mat1->get_diagonal( n, diag );
	EXPECT_EQ( n, 5 );
	int expected1[5] = { 3, 1, 3, 3, 2 };
	EXPECT_ARRAY_EQ( 5, diag, expected1 );

	n = 0;
	mat2->get_diagonal( n, diag );
	EXPECT_EQ( n, 5 );
	int expected2[5] = { 1, 2, 1, 0, 0 };
	EXPECT_ARRAY_EQ( 5, diag, expected2 );

	delete [] diag;
}

TEST_F(MatrixTest,GetOffdiagonalGetsCorrectVectorsForHorizontalRectangularMatrix) {
	vector<int>
			diag   = mat1->get_offdiagonal(0),
			ldiag1 = mat1->get_offdiagonal(-1),
			ldiag2 = mat1->get_offdiagonal(-2),
			ldiag3 = mat1->get_offdiagonal(-3),
			udiag1 = mat1->get_offdiagonal(1),
			udiag2 = mat1->get_offdiagonal(2),
			udiag3 = mat1->get_offdiagonal(3),
			udiag4 = mat1->get_offdiagonal(4),
			udiag5 = mat1->get_offdiagonal(5);
	ASSERT_EQ(diag.size(),5);
	ASSERT_EQ(ldiag1.size(),4);
	ASSERT_EQ(ldiag2.size(),3);
	ASSERT_EQ(ldiag3.size(),2);
	ASSERT_EQ(udiag1.size(),5);
	ASSERT_EQ(udiag2.size(),5);
	ASSERT_EQ(udiag3.size(),5);
	ASSERT_EQ(udiag4.size(),4);
	ASSERT_EQ(udiag5.size(),3);

	EXPECT_ARRAY_EQ( 5, diag,   vector<int>( { 3, 1, 3, 3, 2 }) );
	EXPECT_ARRAY_EQ( 4, ldiag1, vector<int>( { 3, 2, 1, 3 }) );
	EXPECT_ARRAY_EQ( 3, ldiag2, vector<int>( { 0, 3, 3 }) );
	EXPECT_ARRAY_EQ( 2, ldiag3, vector<int>( { 3, 3 }) );
	EXPECT_ARRAY_EQ( 5, udiag1, vector<int>( { 0, 3, 3, 3, 0 }) );
	EXPECT_ARRAY_EQ( 5, udiag2, vector<int>( { 0, 1, 3, 2, 0 }) );
	EXPECT_ARRAY_EQ( 5, udiag3, vector<int>( { 0, 0, 1, 0, 0 }) );
	EXPECT_ARRAY_EQ( 4, udiag4, vector<int>( { 2, 2, 1, 3 }) );
	EXPECT_ARRAY_EQ( 3, udiag5, vector<int>( { 3, 0, 1 }) );
}

TEST_F(MatrixTest,GetOffdiagonalGetsCorrectVectorsForVerticalRectangularMatrix) {
	vector<int>
			diag   = mat2->get_offdiagonal(0),
			ldiag1 = mat2->get_offdiagonal(-1),
			ldiag2 = mat2->get_offdiagonal(-2),
			ldiag3 = mat2->get_offdiagonal(-3),
			ldiag4 = mat2->get_offdiagonal(-4),
			udiag1 = mat2->get_offdiagonal(1),
			udiag2 = mat2->get_offdiagonal(2),
			udiag3 = mat2->get_offdiagonal(3),
			udiag4 = mat2->get_offdiagonal(4);
	ASSERT_EQ(diag.size(),5);
	ASSERT_EQ(ldiag1.size(),5);
	ASSERT_EQ(ldiag2.size(),5);
	ASSERT_EQ(ldiag3.size(),5);
	ASSERT_EQ(udiag1.size(),4);
	ASSERT_EQ(udiag2.size(),3);
	ASSERT_EQ(udiag3.size(),2);
	ASSERT_EQ(udiag4.size(),1);

	EXPECT_ARRAY_EQ( 5, diag,   vector<int>( { 1,2,1,0,0 }) );
	EXPECT_ARRAY_EQ( 5, ldiag1, vector<int>( { 1,0,0,0,0 }) );
	EXPECT_ARRAY_EQ( 5, ldiag2, vector<int>( { 2,2,2,0,1 }) );
	EXPECT_ARRAY_EQ( 5, ldiag3, vector<int>( { 2,1,0,2,1 }) );
	EXPECT_ARRAY_EQ( 4, ldiag4, vector<int>( { 0,0,1,0 }) );
	EXPECT_ARRAY_EQ( 4, udiag1, vector<int>( { 2,1,1,1 }) );
	EXPECT_ARRAY_EQ( 3, udiag2, vector<int>( { 2,2,2 }) );
	EXPECT_ARRAY_EQ( 2, udiag3, vector<int>( { 2,0 }) );
	EXPECT_ARRAY_EQ( 1, udiag4, vector<int>( { 2 }) );
}

TEST_F(MatrixTest,GetOffDiagonalThrowsOutOfRangeForOffsetTooLarge) {
	EXPECT_THROW( {vector<int> v = mat2->get_offdiagonal(5); }, out_of_range );
}

TEST_F(MatrixTest,GetOffdiagonalGetsCorrectArraysForHorizontalRectangularMatrix) {
	int *diag = NCPA::zeros<int>(5);
	size_t n = 0;

	mat1->get_offdiagonal(0,n,diag);
	ASSERT_EQ(n,5);
	EXPECT_ARRAY_EQ( n, diag,   vector<int>( { 3, 1, 3, 3, 2 }) );
	n = 0;

	mat1->get_offdiagonal(-1,n,diag);
	ASSERT_EQ(n,4);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 3, 2, 1, 3 }) );
	n = 0;

	mat1->get_offdiagonal(-2,n,diag);
	ASSERT_EQ(n,3);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 0, 3, 3 }) );
	n = 0;

	mat1->get_offdiagonal(-3,n,diag);
	ASSERT_EQ(n,2);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 3, 3 }) );
	n = 0;

	mat1->get_offdiagonal(1,n,diag);
	ASSERT_EQ(n,5);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 0, 3, 3, 3, 0 } ) );
	n = 0;

	mat1->get_offdiagonal(2,n,diag);
	ASSERT_EQ(n,5);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 0, 1, 3, 2, 0 } ) );
	n = 0;

	mat1->get_offdiagonal(3,n,diag);
	ASSERT_EQ(n,5);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 0, 0, 1, 0, 0 } ) );
	n = 0;

	mat1->get_offdiagonal(4,n,diag);
	ASSERT_EQ(n,4);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 2, 2, 1, 3 } ) );
	n = 0;

	mat1->get_offdiagonal(5,n,diag);
	ASSERT_EQ(n,3);
	EXPECT_ARRAY_EQ( n, diag, vector<int>( { 3, 0, 1 } ) );
	n = 0;
}

TEST_F(MatrixTest,GetOffdiagonalThrowsOutOfRangeForWrongSizeArray) {
	int *v = NCPA::zeros<int>(5);
	size_t n = 5;
	EXPECT_THROW( {mat1->get_offdiagonal(5,n,v); }, out_of_range );
	EXPECT_THROW( {mat1->get_offdiagonal(-1,n,v); }, out_of_range );
	delete [] v;
}

TEST_F(MatrixTest,SetDiagonalSetsDiagonalCorrectly) {
	int v[5] = {1,1,1,1,1};
	mat1->set_diagonal( 5, v );
	EXPECT_ARRAY_EQ( 5, mat1->get_diagonal(), v );
	vector<int> vv( v, v+5 );
	mat2->set_diagonal( vv );
	EXPECT_ARRAY_EQ( 5, mat2->get_diagonal(), v );
	nsmat1->set_diagonal( 1 );
	EXPECT_ARRAY_EQ( 5, nsmat1->get_diagonal(), v );
}

TEST_F(MatrixTest,SetOffdiagonalSetsValuesCorrectly) {
	int v[5] = {1,1,1,1,1};
	mat1->set_offdiagonal( -1, 4, v );
	EXPECT_ARRAY_EQ( 4, mat1->get_offdiagonal(-1), v );
	mat1->set_offdiagonal( 1, 5, v );
	EXPECT_ARRAY_EQ( 5, mat1->get_offdiagonal(1), v );

	vector<int> vv( 4, 1 );
	mat2->set_offdiagonal( 1, vv );
	EXPECT_ARRAY_EQ( 4, mat2->get_offdiagonal(1), vv );
	vv.push_back( 1 );
	mat2->set_offdiagonal( -2, vv );
	EXPECT_ARRAY_EQ( 5, mat2->get_offdiagonal(-2), vv );

	nsmat1->set_offdiagonal( -1, 1 );
	EXPECT_ARRAY_EQ( 4, nsmat1->get_offdiagonal(-1), v );
	nsmat1->set_offdiagonal( 1, 1 );
	EXPECT_ARRAY_EQ( 5, nsmat1->get_offdiagonal(1), v );
}

TEST_F(MatrixTest,SetOffdiagonalThrowsOutOfRangeForWrongSizeArray) {
	int v[5] = {1,1,1,1,1};
	EXPECT_THROW( { mat1->set_offdiagonal( -1, 5, v );}, out_of_range );
	vector<int> vv( v, v+5 );
	EXPECT_THROW( { mat2->set_offdiagonal( 1, vv ); }, out_of_range );
}

