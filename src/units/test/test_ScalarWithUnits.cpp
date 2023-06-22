#include "NCPAUnits.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>

using namespace std;
using namespace NCPA;
using namespace testing;

// test fixture
class ScalarWithUnitsTest : public ::testing::Test {
 protected:

	// Initializations and other setup
	void SetUp() override {
		freezing = ScalarWithUnits( 0.0, CELSIUS );
		boiling  = ScalarWithUnits( 100.0, CELSIUS );
		halfway  = ScalarWithUnits( 50.0, CELSIUS );
		s1 = ScalarWithUnits( boiling );
		s2 = ScalarWithUnits( freezing );
		firstLeg = ScalarWithUnits( 10.0, METERS );
		secondLeg = ScalarWithUnits( 20.0, METERS );
		thirdLeg = ScalarWithUnits( 0.03, KILOMETERS );
	}

	// If there's any cleanup other than normal destructors
	//	void TearDown() override {}

	// class members here
	ScalarWithUnits freezing, boiling, halfway, s0, s1, s2;
	ScalarWithUnits firstLeg, secondLeg, thirdLeg;

	const units_t CELSIUS = Units::fromString("C"),
				  FAHRENHEIT = Units::fromString("F"),
				  KELVIN = Units::fromString("K"),
				  METERS = Units::fromString("m"),
				  KILOMETERS = Units::fromString("km");

};

TEST_F(ScalarWithUnitsTest,ConstructorTests) {
	ASSERT_EQ( s0.get_units(), Units::null() );
	ASSERT_DOUBLE_EQ( s0.get(), 0.0 );

	ASSERT_EQ( boiling.get_units(), CELSIUS );
	ASSERT_DOUBLE_EQ( boiling.get(), 100.0 );

	ASSERT_EQ( s1.get_units(), CELSIUS );
	ASSERT_DOUBLE_EQ( s1.get(), 100.0 );

	// copy constructor
	ScalarWithUnits s4( s0 );
	ASSERT_EQ( s4.get_units(), Units::null() );
	ASSERT_DOUBLE_EQ( s4.get(), 0.0 );
	s0.set( 50.0, CELSIUS );
	ASSERT_EQ( s0.get_units(), CELSIUS );
	ASSERT_DOUBLE_EQ( s0.get(), 50.0 );
	ASSERT_EQ( s4.get_units(), Units::null() );
	ASSERT_DOUBLE_EQ( s4.get(), 0.0 );

	// move constructor?
	s2 = std::move(boiling);
	ASSERT_EQ( s2.get_units(), CELSIUS );
	ASSERT_DOUBLE_EQ( s2.get(), 100.0 );
}

TEST_F(ScalarWithUnitsTest, GetTests) {
	ASSERT_DOUBLE_EQ( freezing.get(), 0.0 );
	ASSERT_EQ( freezing.get_units(), CELSIUS);
}

TEST_F(ScalarWithUnitsTest,GetAsReturnsCorrectValue) {
	ASSERT_DOUBLE_EQ( freezing.get_as(FAHRENHEIT), 32.0 );
	ASSERT_DOUBLE_EQ( freezing.as("Kelvin"), 273.15 );
}

TEST_F(ScalarWithUnitsTest,GetAsThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {double d = freezing.as(METERS);},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, SetValueSetsCorrectValue) {
	freezing.set_value( 212.0 );
	ASSERT_DOUBLE_EQ( freezing.get(), 212.0 );
}

TEST_F(ScalarWithUnitsTest, SetValueDoesNotChangeUnits) {
	units_t u = freezing.get_units();
	freezing.set_value( 212.0 );
	ASSERT_EQ( freezing.get_units(), u );
}

TEST_F(ScalarWithUnitsTest, SetUnitsSetsCorrectUnits) {
	freezing.set_units( FAHRENHEIT );
	ASSERT_EQ( freezing.get_units(), FAHRENHEIT );
}

TEST_F(ScalarWithUnitsTest, SetUnitsDoesNotChangeValue) {
	freezing.set_units( FAHRENHEIT );
	ASSERT_DOUBLE_EQ( freezing.get(), 0.0 );
}

TEST_F(ScalarWithUnitsTest, SetUnitsWorksWithString) {
	freezing.set_units( "F" );
	ASSERT_EQ( freezing.get_units(), FAHRENHEIT );
}

TEST_F(ScalarWithUnitsTest, SetOverwritesBothValueAndUnits) {
	// Set units to 32 degrees fahrenheit
	firstLeg.set( 32.0, FAHRENHEIT );
	ASSERT_DOUBLE_EQ( firstLeg.get(), 32.0 );
	ASSERT_EQ( firstLeg.get_units(), FAHRENHEIT );
}

TEST_F(ScalarWithUnitsTest, ConvertUnitsConvertsCorrectly) {
	freezing.convert_units(FAHRENHEIT);
	ASSERT_DOUBLE_EQ(freezing.get(),32.0);
}

TEST_F(ScalarWithUnitsTest,ConvertUnitsThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {freezing.convert_units(METERS);},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, EqualityOperatorIsTrueForEqual) {
	ASSERT_TRUE( s1 == boiling );
}

TEST_F(ScalarWithUnitsTest, EqualityOperatorIsFalseForUnequal) {
	ASSERT_FALSE( s1 == freezing );
}

TEST_F(ScalarWithUnitsTest, EqualityOperatorIsTrueForEqualWithNullUnits) {
	boiling.set_units( NCPA::Units::null() );
	ASSERT_TRUE( s1 == boiling );
}

TEST_F(ScalarWithUnitsTest, EqualityOperatorWorksWithDouble) {
	ASSERT_TRUE( s1 == 100.0 );
}

TEST_F(ScalarWithUnitsTest, EqualityOperatorCommutesWithDouble) {
	ASSERT_TRUE( 100.0 == s1 );
}

TEST_F(ScalarWithUnitsTest, InequalityOperatorIsTrueForUnequal) {
	ASSERT_TRUE( freezing != boiling );
}

TEST_F(ScalarWithUnitsTest, InequalityOperatorIsFalseForEqual) {
	ASSERT_FALSE( boiling != boiling );
}

TEST_F(ScalarWithUnitsTest, InequalityOperatorIsFalseForEqualWithNullUnits) {
	boiling.set_units( Units::null() );
	ASSERT_FALSE( freezing == boiling );
}

TEST_F(ScalarWithUnitsTest, InequalityOperatorWorksWithDouble) {
	ASSERT_TRUE( boiling != 0.0 );
}

TEST_F(ScalarWithUnitsTest, InequalityOperatorCommutesWithDouble) {
	ASSERT_TRUE( 0.0 != boiling );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOperatorTrueForGreater) {
	ASSERT_TRUE( boiling > freezing );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOperatorFalseForLess) {
	ASSERT_FALSE( freezing > boiling );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOperatorFalseForEqual) {
	ASSERT_FALSE( boiling > boiling );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOperatorTrueForLesserDouble) {
	ASSERT_TRUE( boiling > 0.0 );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOperatorFalseForLesserDouble) {
	ASSERT_FALSE( boiling > 200.0 );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOperatorFalseForEqualDouble) {
	ASSERT_FALSE( freezing > 0.0 );
}

TEST_F(ScalarWithUnitsTest,GreaterThanOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {bool b = freezing > firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, LessThanOperatorFalseForGreater) {
	ASSERT_FALSE( boiling < freezing );
}

TEST_F(ScalarWithUnitsTest, LessThanOperatorTrueForLesser) {
	ASSERT_TRUE( freezing < boiling );
}

TEST_F(ScalarWithUnitsTest, LessThanOperatorFalseForEqual) {
	ASSERT_FALSE( boiling < boiling );
}

TEST_F(ScalarWithUnitsTest, LessThanOperatorTrueForGreaterDouble) {
	ASSERT_TRUE( boiling < 200.0 );
}

TEST_F(ScalarWithUnitsTest, LessThanOperatorFalseForLesserDouble) {
	ASSERT_FALSE( boiling < 0.0 );
}

TEST_F(ScalarWithUnitsTest, LessThanOperatorFalseForEqualDouble) {
	ASSERT_FALSE( freezing < 0.0 );
}

TEST_F(ScalarWithUnitsTest,LessThanOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {bool b = freezing < firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOrEqualOperatorTrueForGreater) {
	ASSERT_TRUE( boiling >= freezing );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOrEqualOperatorFalseForLess) {
	ASSERT_FALSE( freezing >= boiling );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOrEqualOperatorTrueForEqual) {
	ASSERT_TRUE( boiling >= boiling );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOrEqualOperatorTrueForLesserDouble) {
	ASSERT_TRUE( boiling >= 0.0 );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOrEqualOperatorFalseForLesserDouble) {
	ASSERT_FALSE( boiling >= 200.0 );
}

TEST_F(ScalarWithUnitsTest, GreaterThanOrEqualOperatorTrueForEqualDouble) {
	ASSERT_TRUE( freezing >= 0.0 );
}

TEST_F(ScalarWithUnitsTest,GreaterThanOrEqualOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {bool b = freezing >= firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, LessThanOrEqualOperatorFalseForGreater) {
	ASSERT_FALSE( boiling <= freezing );
}

TEST_F(ScalarWithUnitsTest, LessThanOrEqualOperatorTrueForLesser) {
	ASSERT_TRUE( freezing <= boiling );
}

TEST_F(ScalarWithUnitsTest, LessThanOrEqualOperatorTrueForEqual) {
	ASSERT_TRUE( boiling <= boiling );
}

TEST_F(ScalarWithUnitsTest, LessThanOrEqualOperatorTrueForGreaterDouble) {
	ASSERT_TRUE( boiling <= 200.0 );
}

TEST_F(ScalarWithUnitsTest, LessThanOrEqualOperatorFalseForLesserDouble) {
	ASSERT_FALSE( boiling <= 0.0 );
}

TEST_F(ScalarWithUnitsTest, LessThanOrEqualOperatorTrueForEqualDouble) {
	ASSERT_TRUE( freezing <= 0.0 );
}

TEST_F(ScalarWithUnitsTest,LessThanOrEqualOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {bool b = freezing <= firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, NegationOperatorMakesValueNegative) {
	ASSERT_DOUBLE_EQ( (-boiling).get(), -(boiling.get()) );
}

TEST_F(ScalarWithUnitsTest, NegationOperatorDoesNotChangeUnits) {
	ASSERT_EQ( (-boiling).get_units(), boiling.get_units() );
}

TEST_F(ScalarWithUnitsTest, PlusOperatorGivesCorrectSum) {
	ScalarWithUnits firstAndSecond = firstLeg + secondLeg;
	ASSERT_DOUBLE_EQ( firstAndSecond.get(), 30.0 );
}

TEST_F(ScalarWithUnitsTest, PlusOperatorGivesCorrectSumWithUnitConversion) {
	ScalarWithUnits secondAndThird = secondLeg + thirdLeg;
	ASSERT_DOUBLE_EQ( secondAndThird.get(), 50.0 );
}

TEST_F(ScalarWithUnitsTest, PlusOperatorKeepsCorrectUnitsWithUnitConversion) {
	ScalarWithUnits secondAndThird = secondLeg + thirdLeg;
	ASSERT_EQ( secondAndThird.get_units(), METERS );
}

TEST_F(ScalarWithUnitsTest, PlusOperatorChainsCorrectly) {
	ScalarWithUnits fullTrip = firstLeg + secondLeg + thirdLeg;
	ASSERT_DOUBLE_EQ( fullTrip.get(), 60.0 );
}

TEST_F(ScalarWithUnitsTest, PlusOperatorChainKeepsCorrectUnits) {
	ScalarWithUnits fullTrip = firstLeg + secondLeg + thirdLeg;
	ASSERT_EQ( fullTrip.get_units(), METERS );
}

TEST_F(ScalarWithUnitsTest, PlusOperatorWorksWithDouble) {
	ScalarWithUnits superHeated = boiling + 50.0;
	ASSERT_DOUBLE_EQ( superHeated.get(), 150.0 );
	ASSERT_EQ( superHeated.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest,PlusOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {ScalarWithUnits b = freezing + firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, MinusOperatorGivesCorrectDifference) {
	ScalarWithUnits firstAndSecond = firstLeg - secondLeg;
	ASSERT_DOUBLE_EQ( firstAndSecond.get(), -10.0 );
}

TEST_F(ScalarWithUnitsTest, MinusOperatorGivesCorrectDifferenceWithUnitConversion) {
	ScalarWithUnits secondAndThird = secondLeg - thirdLeg;
	ASSERT_DOUBLE_EQ( secondAndThird.get(), -10.0 );
}

TEST_F(ScalarWithUnitsTest, MinusOperatorKeepsCorrectUnitsWithUnitConversion) {
	ScalarWithUnits secondAndThird = secondLeg - thirdLeg;
	ASSERT_EQ( secondAndThird.get_units(), METERS );
}

TEST_F(ScalarWithUnitsTest, MinusOperatorChainsCorrectly) {
	ScalarWithUnits fullTrip = firstLeg - secondLeg - thirdLeg;
	ASSERT_DOUBLE_EQ( fullTrip.get(), -40.0 );
}

TEST_F(ScalarWithUnitsTest, MinusOperatorChainKeepsCorrectUnits) {
	ScalarWithUnits fullTrip = firstLeg - secondLeg - thirdLeg;
	ASSERT_EQ( fullTrip.get_units(), METERS );
}

TEST_F(ScalarWithUnitsTest, MinusOperatorWorksWithDouble) {
	ScalarWithUnits superCooled = freezing - 50.0;
	ASSERT_DOUBLE_EQ( superCooled.get(), -50.0 );
	ASSERT_EQ( superCooled.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest,MinusOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {ScalarWithUnits b = freezing - firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, PlusEqualsOperatorGivesCorrectSum) {
	firstLeg += secondLeg;
	ASSERT_DOUBLE_EQ( firstLeg.get(), 30.0 );
}

TEST_F(ScalarWithUnitsTest, PlusEqualsOperatorGivesCorrectSumWithUnitConversion) {
	secondLeg += thirdLeg;
	ASSERT_DOUBLE_EQ( secondLeg.get(), 50.0 );
}

TEST_F(ScalarWithUnitsTest, PlusEqualsOperatorKeepsCorrectUnitsWithUnitConversion) {
	secondLeg += thirdLeg;
	ASSERT_EQ( secondLeg.get_units(), METERS );
}


TEST_F(ScalarWithUnitsTest, PlusEqualsOperatorWorksWithDouble) {
	boiling += 50.0;
	ASSERT_DOUBLE_EQ( boiling.get(), 150.0 );
	ASSERT_EQ( boiling.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest,PlusEqualsOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {freezing += firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, MinusEqualsOperatorGivesCorrectDifference) {
	firstLeg -= secondLeg;
	ASSERT_DOUBLE_EQ( firstLeg.get(), -10.0 );
}

TEST_F(ScalarWithUnitsTest, MinusEqualsOperatorGivesCorrectDifferenceWithUnitConversion) {
	secondLeg -= thirdLeg;
	ASSERT_DOUBLE_EQ( secondLeg.get(), -10.0 );
}

TEST_F(ScalarWithUnitsTest, MinusEqualsOperatorWorksWithDouble) {
	freezing -= 50.0;
	ASSERT_DOUBLE_EQ( freezing.get(), -50.0 );
	ASSERT_EQ( freezing.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest,MinusEqualsOperatorThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {freezing -= firstLeg;},
				invalid_conversion );
}

TEST_F(ScalarWithUnitsTest, MultiplyEqualsOperatorScalesCorrectly) {
	boiling *= 3.0;
	ASSERT_DOUBLE_EQ( boiling.get(), 300.0 );
	ASSERT_EQ( boiling.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest, DivideEqualsOperatorScalesCorrectly) {
	boiling /= 4.0;
	ASSERT_DOUBLE_EQ( boiling.get(), 25.0 );
	ASSERT_EQ( boiling.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest, MultiplyOperatorScalesCorrectly) {
	ScalarWithUnits temp = boiling * 3.0;
	ASSERT_DOUBLE_EQ( temp.get(), 300.0 );
	ASSERT_EQ( temp.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest, DivideOperatorScalesCorrectly) {
	ScalarWithUnits temp = boiling / 4.0;
	ASSERT_DOUBLE_EQ( temp.get(), 25.0 );
	ASSERT_EQ( temp.get_units(), CELSIUS );
}

TEST_F(ScalarWithUnitsTest, OverGivesCorrectRatio) {
	ASSERT_DOUBLE_EQ( boiling.over(halfway), 2.0 );
}

TEST_F(ScalarWithUnitsTest,OverThrowsOutOfRangeOnInvalidConversion) {
	EXPECT_THROW( {double d = boiling.over(firstLeg);},
				invalid_conversion );
}



//EXPECT_THROW( {out = Units::convert( in, units_t::PRESSURE_MILLIBARS, units_t::ANGLE_DEGREES );},
//			invalid_conversion );
