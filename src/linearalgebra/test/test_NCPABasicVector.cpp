#include "NCPACommon.h"
#include "NCPAVector.h"
#include "NCPABasicVector.h"
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
class NCPABasicVectorTest : public ::testing::Test {
protected:

	// Initializations and other setup
	void SetUp() override {
		inds_a = new size_t[4];
		vals_a = new double[4];
		for (size_t i = 0; i < 4; i++) {
			inds_a[i] = i;
			vals_a[i] = (double)(i+3);
		}
		inds_off_a = new size_t[2];
		inds_off_a[0] = 1;
		inds_off_a[1] = 2;

		a = new NCPA::BasicVector<double>();
		b = new NCPA::BasicVector<double>( 4 );
		c = new NCPA::BasicVector<double>( 4, inds, avals );
		d = new NCPA::BasicVector<double>( 4, 4, inds_a, vals_a );
		e = new NCPA::BasicVector<double>( 4, { 0, 1, 2, 3 }, { 3, 4, 5, 6 } );
		f = new NCPA::BasicVector<double>( 4, 5.0 );
		g = new NCPA::BasicVector<double>( 4, inds_off, 5.0 );
		h = new NCPA::BasicVector<double>( 4, 2, inds_off_a, 5.0 );
		m = new NCPA::BasicVector<double>( 4, { 1, 2 }, 5.0 );
		two = new NCPA::BasicVector<double>( 2 );
		q = new NCPA::BasicVector<double>( 4, inds_off, avals );
		r = new NCPA::BasicVector<double>( 4, inds, bvals );
	}

	// If there's any cleanup other than normal destructors
	//void TearDown() override {}

	// class members here
	NCPA::Vector<double> *a, *b, *c, *d, *e, *f, *g, *h, *m, *two, *q, *r;
	std::vector<size_t> inds = { 0, 1, 2, 3 };
	std::vector<size_t> inds_off = { 1, 2 };
	std::vector<double> avals = { 3, 4, 5, 6 };
	std::vector<double> bvals = { -2, -3, 4, 5 };
	size_t *inds_a, *inds_off_a;
	double *vals_a;
};

TEST_F(NCPABasicVectorTest,BuildWithNoArgumentsBuildsZeroLengthVector) {
	EXPECT_EQ( a->size(), 0 );
}

TEST_F(NCPABasicVectorTest,BuildWithNoElementsBuildsEmptyVector) {
	EXPECT_EQ( b->size(), 4 );
	for (size_t i = 0; i < b->size(); i++) {
		EXPECT_DOUBLE_EQ( b->get(i), 0.0 );
	}
}

TEST_F(NCPABasicVectorTest,BuildWithStdVectorsBuildsCorrectly) {
	EXPECT_EQ( c->size(), 4 );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_DOUBLE_EQ( c->get(i), avals[i] );
	}
}

TEST_F(NCPABasicVectorTest,BuildWithArraysBuildsCorrectly) {
	EXPECT_EQ( d->size(), 4 );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_DOUBLE_EQ( d->get(i), vals_a[i] );
	}
}

TEST_F(NCPABasicVectorTest,BuildWithInitializerListsBuildsCorrectly) {
	EXPECT_EQ( e->size(), 4 );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_DOUBLE_EQ( e->get(i), vals_a[i] );
	}
}

TEST_F(NCPABasicVectorTest,BuildWithConstantBuildsCorrectly) {
	EXPECT_EQ( f->size(), 4 );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_DOUBLE_EQ( f->get(i), 5.0 );
	}
}

TEST_F(NCPABasicVectorTest,BuildWithConstantAndIndicesBuildsCorrectly) {
	EXPECT_EQ( g->size(), 4 );
	EXPECT_DOUBLE_EQ( g->get(0), 0.0 );
	EXPECT_DOUBLE_EQ( g->get(1), 5.0 );
	EXPECT_DOUBLE_EQ( g->get(2), 5.0 );
	EXPECT_DOUBLE_EQ( g->get(3), 0.0 );
}

TEST_F(NCPABasicVectorTest,BuildWithConstantAndArrayOfIndicesBuildsCorrectly) {
	EXPECT_EQ( h->size(), 4 );
	EXPECT_DOUBLE_EQ( h->get(0), 0.0 );
	EXPECT_DOUBLE_EQ( h->get(1), 5.0 );
	EXPECT_DOUBLE_EQ( h->get(2), 5.0 );
	EXPECT_DOUBLE_EQ( h->get(3), 0.0 );
}

TEST_F(NCPABasicVectorTest,BuildWithConstantAndListOfIndicesBuildsCorrectly) {
	EXPECT_EQ( m->size(), 4 );
	EXPECT_DOUBLE_EQ( m->get(0), 0.0 );
	EXPECT_DOUBLE_EQ( m->get(1), 5.0 );
	EXPECT_DOUBLE_EQ( m->get(2), 5.0 );
	EXPECT_DOUBLE_EQ( m->get(3), 0.0 );
}


TEST_F(NCPABasicVectorTest,BuildWithStdVectorsAndMissingIndicesBuildsCorrectly) {
	EXPECT_EQ( q->size(), 4 );
	EXPECT_DOUBLE_EQ( q->get( 0 ), 0.0 );
	EXPECT_DOUBLE_EQ( q->get( 1 ), 3.0 );
	EXPECT_DOUBLE_EQ( q->get( 2 ), 4.0 );
	EXPECT_DOUBLE_EQ( q->get( 3 ), 0.0 );
}

TEST_F(NCPABasicVectorTest,ClearClearsVector) {
	m->clear();
	EXPECT_EQ(m->size(), 0);
	EXPECT_THROW( { m->get(0); }, out_of_range );
}

TEST_F(NCPABasicVectorTest,ClearClearsElement) {
	q->clear( 1 );
	EXPECT_EQ( q->size(), 4 );
	EXPECT_DOUBLE_EQ( q->get( 1 ), 0.0 );
}

TEST_F(NCPABasicVectorTest,ResizeShrinksVector) {
	m->resize(2);
	EXPECT_EQ(m->size(),2);
	EXPECT_DOUBLE_EQ( m->get(0), 0.0 );
	EXPECT_DOUBLE_EQ( m->get(1), 5.0 );
	EXPECT_THROW( { m->get(2); }, out_of_range );
}

TEST_F(NCPABasicVectorTest,ResizeGrowsVector) {
	m->resize( 6 );
	EXPECT_EQ(m->size(),6);
	EXPECT_DOUBLE_EQ( m->get(5), 0.0 );
}

TEST_F(NCPABasicVectorTest,CloneWorks) {
	NCPA::Vector<double> *temp = c->clone();
	EXPECT_EQ( c->size(), temp->size() );
	for (size_t i = 0; i < c->size(); i++) {
		EXPECT_DOUBLE_EQ( c->get(i), temp->get(i) );
	}
	size_t n = c->size();
	delete c;
	EXPECT_EQ( temp->size(), n );
}

TEST_F(NCPABasicVectorTest,SetWorks) {
	m->set(0,13.0);
	EXPECT_DOUBLE_EQ( m->get(0), 13.0 );
}

TEST_F(NCPABasicVectorTest,ScaleWorksCorrectly) {
	q->scale( -2.0 );
	EXPECT_DOUBLE_EQ( q->get( 0 ), 0.0 );
	EXPECT_DOUBLE_EQ( q->get( 1 ), -6.0 );
	EXPECT_DOUBLE_EQ( q->get( 2 ), -8.0 );
	EXPECT_DOUBLE_EQ( q->get( 3 ), 0.0 );
}

TEST_F(NCPABasicVectorTest,AddWorksCorrectly) {
	q->add( c );
	EXPECT_DOUBLE_EQ( q->get( 0 ), 3.0 );
	EXPECT_DOUBLE_EQ( q->get( 1 ), 7.0 );
	EXPECT_DOUBLE_EQ( q->get( 2 ), 9.0 );
	EXPECT_DOUBLE_EQ( q->get( 3 ), 6.0 );
}

TEST_F(NCPABasicVectorTest,SubtractWorksCorrectly) {
	q->subtract(c);
	EXPECT_DOUBLE_EQ( q->get( 0 ), -3.0 );
	EXPECT_DOUBLE_EQ( q->get( 1 ), -1.0 );
	EXPECT_DOUBLE_EQ( q->get( 2 ), -1.0 );
	EXPECT_DOUBLE_EQ( q->get( 3 ), -6.0 );
}

TEST_F(NCPABasicVectorTest,ScalarProductWorksCorrectly) {
	double prod = c->scalar_product( r );
	EXPECT_DOUBLE_EQ( prod, 32.0 );
	prod = r->scalar_product( c );
	EXPECT_DOUBLE_EQ( prod, 32.0 );
}

TEST_F(NCPABasicVectorTest,ScalarProductThrowsOutOfRangeIfSizesDoNotMatch) {
	EXPECT_THROW( { c->scalar_product( two ); }, out_of_range );
}

TEST_F(NCPABasicVectorTest,ProductMultipliesCorrectly) {
	Vector<double> *prod = c->product( r );
	for (size_t i = 0; i < r->size(); i++) {
		EXPECT_DOUBLE_EQ( prod->get(i), c->get(i) * r->get(i) );
	}
}

TEST_F(NCPABasicVectorTest,AsArrayReturnsCorrectly) {
	double *buffer = new double[ 4 ];
	c->as_array( 4, buffer );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_DOUBLE_EQ( c->get(i), buffer[i] );
	}
	delete [] buffer;
}

TEST_F(NCPABasicVectorTest,AsArrayAllocatesMemory) {
	double *buffer = nullptr;
	c->as_array( 4, buffer );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_DOUBLE_EQ( c->get(i), buffer[i] );
	}
	delete [] buffer;
}

TEST_F(NCPABasicVectorTest,AsStdReturnsCorrectly) {
	EXPECT_EQ( q->as_std(), std::vector<double>( { 0, 3, 4, 0 } ) );
}

TEST_F(NCPABasicVectorTest,FromStdBuildsCorrectVector) {
	two->from_std( avals );
	EXPECT_EQ( two->size(), 4 );
	for (size_t i = 0; i < 4; i++) {
//		cout << two->get(i) << endl;
		EXPECT_DOUBLE_EQ( two->get(i), avals[i] );
	}
}

TEST_F(NCPABasicVectorTest,GetDefinedIndicesWorksCorrectly) {
	std::vector<size_t> i1, i3;
	i1 = c->get_defined_indices();
	EXPECT_EQ( i1.size(), 4 );
	for (size_t i = 0; i < 4; i++) {
		EXPECT_EQ( i1[i], i );
	}
	i3 = g->get_defined_indices();
	EXPECT_EQ( i3.size(), 2 );
	for (size_t i = 0; i < 2; i++) {
		EXPECT_EQ( i3[i], i+1 );
	}
}

// Tests that don't use the fixture
//TEST(NCPABasicVectorTest,TestName2) {
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
