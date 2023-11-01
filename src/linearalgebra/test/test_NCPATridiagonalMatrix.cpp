#include "NCPALinearAlgebra.h"
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <complex>

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
class TridiagonalMatrixTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		tmat1 = new TridiagonalMatrix<int>( 7, 7 );
		tmat1->set_diagonal(1);
		tmat1->set_offdiagonal(1,2);
		tmat1->set_offdiagonal(-1,3);

		tmat2 = new TridiagonalMatrix<int>( *tmat1 );

		smat2 = new SparseMatrix<int>( 7, 7 );
		smat2->set_diagonal(13);
		smat2->set_offdiagonal( -1, 6 );
		smat2->set_offdiagonal( -2, 9 );
		smat2->set_offdiagonal( 1, 4 );
		smat2->set_offdiagonal( 2, 4 );
		smat2->set( 0, 0, 7 );
		smat2->set( 6, 6, 7 );

		mat2 = static_cast<Matrix<int>*>(smat2);

		dmat = new TridiagonalMatrix<double>( 7, 7 );
		dmat->set_diagonal(1.0)->set_offdiagonal(1,2.0)->set_offdiagonal(-1,3.0);

		complex<double> cd( 1.0, -2.0 ), ud( 2.0, 0.5 ), ld( 3.0, 0.0 );
		cmat = new TridiagonalMatrix<complex<double>>( 7, 7 );
		cmat->set_diagonal( cd )->set_offdiagonal(1,ud)->set_offdiagonal(-1,ld);

	}

	// If there's any cleanup other than normal destructors
	//void TearDown() override {}

	// class members here
	TridiagonalMatrix<int> *tmat1, *tmat2;
	Matrix<int> *mat1, *mat2;
	SparseMatrix<int> *smat2;
	TridiagonalMatrix<double> *dmat;
	TridiagonalMatrix<complex<double>> *cmat;

};

TEST_F(TridiagonalMatrixTest,EqualityOperatorsWorkCorrectly) {
	TridiagonalMatrix<int> *t = new TridiagonalMatrix<int>( 7, 7 );
	vector<int> ldiag( 6, 3 ), udiag( 6, 2 ), diag( 7, 1 );
	t->set_diagonal(1);
	t->set_offdiagonal(1,2);
	t->set_offdiagonal(-1,3);
	EXPECT_TRUE( *tmat1 == *t );
	EXPECT_FALSE( *tmat1 != *t );
	tmat1->set( 2, 3, 0 );
	EXPECT_FALSE( *tmat1 == *t );
	EXPECT_TRUE( *tmat1 != *t );
}

TEST_F(TridiagonalMatrixTest,CopyConstructorCreatesCorrectCopy) {
	EXPECT_TRUE( *tmat2 == *tmat1 );
}

TEST_F(TridiagonalMatrixTest,AssignmentOperatorCreatesCopy) {
	TridiagonalMatrix<int> tmat5 = *tmat1;
	EXPECT_TRUE( tmat5 == *tmat1 );
	EXPECT_FALSE( &tmat5 == tmat1 );
}

TEST_F(TridiagonalMatrixTest,CloneMethodCreatesCopy) {
	Matrix<int> *mat1 = tmat2->clone();
	EXPECT_TRUE( *mat1 == *tmat2 );
	EXPECT_FALSE( mat1 == tmat2 );
	delete mat1;
}

TEST_F(TridiagonalMatrixTest,SetMethodSetsCorrectly) {
	tmat1->set( 1, 2, 5 );
	EXPECT_EQ( tmat1->get( 1, 2 ), 5 );
}

TEST_F(TridiagonalMatrixTest,SetMethodThrowsOutOfRangeIfTooFarOffDiagonal) {
	EXPECT_THROW( { tmat1->set( 1, 3, 5 ); }, out_of_range );
}

TEST_F(TridiagonalMatrixTest,SetRowWorksCorrectlyWithArrays) {
	size_t inds[3] = { 0, 1, 2 };
	int vals[3] = { 5, 6, 5 };
	tmat1->set_row( 1, 3, inds, vals );
	EXPECT_EQ( tmat1->get( 1, 0 ), 5 );
	EXPECT_EQ( tmat1->get( 1, 1 ), 6 );
	EXPECT_EQ( tmat1->get( 1, 2 ), 5 );
	EXPECT_EQ( tmat1->get( 1, 3 ), 0 );
}

TEST_F(TridiagonalMatrixTest,SetRowWorksCorrectlyWithIndicesOutOfOrder) {
	size_t inds[3] = { 1, 2, 0 };
	int vals[3] = { 6, 5, 5 };
	tmat1->set_row( 1, 3, inds, vals );
	EXPECT_EQ( tmat1->get( 1, 0 ), 5 );
	EXPECT_EQ( tmat1->get( 1, 1 ), 6 );
	EXPECT_EQ( tmat1->get( 1, 2 ), 5 );
}

TEST_F(TridiagonalMatrixTest,SetRowThrowsOutOfRangeIfTooFarOffDiagonal) {
	size_t inds[3] = { 2, 3, 4 };
	int vals[3] = { 5, 6, 5 };
	EXPECT_THROW( { tmat1->set_row( 1, 3, inds, vals ); }, out_of_range );
}

TEST_F(TridiagonalMatrixTest,SetRowWorksCorrectlyFor3ElementVector) {
	vector<int> vals( {5, 6, 5} );
	tmat1->set_row( 3, vals );
	EXPECT_EQ( tmat1->get( 3, 2 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 3 ), 6 );
	EXPECT_EQ( tmat1->get( 3, 4 ), 5 );
}

TEST_F(TridiagonalMatrixTest,SetRowWorksCorrectlyForFullLengthVector) {
	vector<int> vals( { 1, 2, 5, 6, 5, 2, 1 } );
	tmat1->set_row( 3, vals );
	EXPECT_EQ( tmat1->get( 3, 1 ), 0 );
	EXPECT_EQ( tmat1->get( 3, 2 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 3 ), 6 );
	EXPECT_EQ( tmat1->get( 3, 4 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 5 ), 0 );
}

TEST_F(TridiagonalMatrixTest,SetRowWorksCorrectlyWithTwoVectors) {
	vector<size_t> inds( { 2, 3, 4 } );
	vector<int> vals( { 5, 6, 5 } );
	tmat1->set_row( 3, inds, vals );
	EXPECT_EQ( tmat1->get( 3, 1 ), 0 );
	EXPECT_EQ( tmat1->get( 3, 2 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 3 ), 6 );
	EXPECT_EQ( tmat1->get( 3, 4 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 5 ), 0 );
}

TEST_F(TridiagonalMatrixTest,SetColumnWorksCorrectlyWithArrays) {
	size_t inds[3] = { 0, 1, 2 };
	int vals[3] = { 5, 6, 5 };
	tmat1->set_column( 1, 3, inds, vals );
	EXPECT_EQ( tmat1->get( 0, 1 ), 5 );
	EXPECT_EQ( tmat1->get( 1, 1 ), 6 );
	EXPECT_EQ( tmat1->get( 2, 1 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 1 ), 0 );
}

TEST_F(TridiagonalMatrixTest,SetColumnWorksCorrectlyWithIndicesOutOfOrder) {
	size_t inds[3] = { 1, 2, 0 };
	int vals[3] = { 6, 5, 5 };
	tmat1->set_column( 1, 3, inds, vals );
	EXPECT_EQ( tmat1->get( 0, 1 ), 5 );
	EXPECT_EQ( tmat1->get( 1, 1 ), 6 );
	EXPECT_EQ( tmat1->get( 2, 1 ), 5 );
}

TEST_F(TridiagonalMatrixTest,SetColumnThrowsOutOfRangeIfTooFarOffDiagonal) {
	size_t inds[3] = { 2, 3, 4 };
	int vals[3] = { 5, 6, 5 };
	EXPECT_THROW( { tmat1->set_column( 1, 3, inds, vals ); }, out_of_range );
}

TEST_F(TridiagonalMatrixTest,SetColumnWorksCorrectlyFor3ElementVector) {
	vector<int> vals( {5, 6, 5} );
	tmat1->set_column( 3, vals );
	EXPECT_EQ( tmat1->get( 2, 3 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 3 ), 6 );
	EXPECT_EQ( tmat1->get( 4, 3 ), 5 );
}

TEST_F(TridiagonalMatrixTest,SetColumnWorksCorrectlyForFullLengthVector) {
	vector<int> vals( { 1, 2, 5, 6, 5, 2, 1 } );
	tmat1->set_column( 3, vals );
	EXPECT_EQ( tmat1->get( 1, 3 ), 0 );
	EXPECT_EQ( tmat1->get( 2, 3 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 3 ), 6 );
	EXPECT_EQ( tmat1->get( 4, 3 ), 5 );
	EXPECT_EQ( tmat1->get( 5, 3 ), 0 );
}

TEST_F(TridiagonalMatrixTest,SetColumnWorksCorrectlyWithTwoVectors) {
	vector<size_t> inds( { 2, 3, 4 } );
	vector<int> vals( { 5, 6, 5 } );
	tmat1->set_column( 3, inds, vals );
	EXPECT_EQ( tmat1->get( 1, 3 ), 0 );
	EXPECT_EQ( tmat1->get( 2, 3 ), 5 );
	EXPECT_EQ( tmat1->get( 3, 3 ), 6 );
	EXPECT_EQ( tmat1->get( 4, 3 ), 5 );
	EXPECT_EQ( tmat1->get( 5, 3 ), 0 );
}

TEST_F(TridiagonalMatrixTest,MultiplyWorksCorrectly) {
	Matrix<int> *product = new SparseMatrix<int>(7,7);
	TridiagonalMatrix<int>::multiply( tmat1, tmat1, product );
	EXPECT_EQ( *product, *mat2 );

	Matrix<int> *product2 = TridiagonalMatrix<int>::multiply( tmat1, tmat1 );
	EXPECT_EQ( *product2, *mat2 );
}

TEST_F(TridiagonalMatrixTest,GetColumnIndicesReturnsCorrectIndices) {
	EXPECT_ARRAY_EQ( 2, tmat1->get_column_indices( 0 ), vector<size_t>( { 0, 1 } ) );
	EXPECT_ARRAY_EQ( 3, tmat1->get_column_indices( 1 ), vector<size_t>( { 0, 1, 2 } ) );
	EXPECT_ARRAY_EQ( 3, tmat1->get_column_indices( 4 ), vector<size_t>( { 3, 4, 5 } ) );
	EXPECT_ARRAY_EQ( 2, tmat1->get_column_indices( 6 ), vector<size_t>( { 5, 6 } ) );
}

TEST_F(TridiagonalMatrixTest,MatrixSystemSolverWorksCorrectly) {
	std::vector<double> y( { 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0 } );
	std::vector<double> expected( {
		0.0538033,
		0.4730983,
		0.6827458,
		0.4489796,
		0.7513915,
		0.4508349,
		-0.3525046
	} );

	std::vector<double> x = dmat->solve( y );
	EXPECT_DOUBLE_ARRAY_NEAR( 7, x, expected, 1e-5 );
//	for (auto it = x.cbegin(); it != x.cend(); ++it) {
//		cout << *it << endl;
//	}
}

TEST_F(TridiagonalMatrixTest,MatrixSystemSolverWorksCorrectlyForComplexValues) {
//	cmat->print();

	complex<double> base( 1.0, 1.0 );
	std::vector<complex<double>> y( {
		1.0*base,
		2.0*base,
		3.0*base,
		4.0*base,
		5.0*base,
		6.0*base,
		7.0*base
	} );
	std::vector<complex<double>> expected( {
		complex<double>( 0.1577105, 0.2277261 ),
		complex<double>( 0.3100052, 0.4663462 ),
		complex<double>( 0.3067259, 0.6585615 ),
		complex<double>( 0.3929875, 0.6796790 ),
		complex<double>( 0.8753552, 0.8464670 ),
		complex<double>( 1.0442588, 1.6715385 ),
		complex<double>(-0.0207091, 1.9439664 ),
	} );

	std::vector<complex<double>> x = cmat->solve( y );

	for (size_t i = 0; i < 7; i++) {
		EXPECT_NEAR( x[i].real(), expected[i].real(), 1e-5 );
		EXPECT_NEAR( x[i].imag(), expected[i].imag(), 1e-5 );
//		cout << "for y = " << y[i] << ", x = " << x[i]
//			<< ", expected = " << expected[i] << endl;
	}
}

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
