#include "NCPACommon.h"
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

#ifndef EXPECT_DOUBLE_ARRAY_NEAR(N,A,Ex,T)
#define EXPECT_DOUBLE_ARRAY_NEAR(N,A,Ex,T) for (size_t i=0; i<N; i++) { EXPECT_NEAR(A[i],Ex[i],T); }
#endif

#ifndef EXPECT_DOUBLE_ARRAY_EQ(N,A,Ex)
#define EXPECT_DOUBLE_ARRAY_EQ(N,A,Ex) for (size_t i=0; i<N; i++) { EXPECT_DOUBLE_EQ(A[i],Ex[i]); }
#endif

class Atmosphere2DTest : public ::testing::Test {
protected:
	// Initializations and other setup
	void SetUp() override {

	}

	// If there's any cleanup other than normal destructors
	//	void TearDown() override {}

	// class members here
	Atmosphere2D a0, a1;
};

TEST_F(Atmosphere2DTest,DefaultConstructorIsEmpty) {

}
