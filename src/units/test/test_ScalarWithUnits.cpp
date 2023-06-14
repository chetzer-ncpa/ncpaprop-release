#include "NCPAUnits.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <utility>

using namespace std;
using namespace NCPA;
using namespace testing;

TEST(ScalarWithUnitsTest,ConstructorTests) {
	ScalarWithUnits s1;
	ASSERT_EQ( s1.get_units(), UNITS_NONE );
	ASSERT_DOUBLE_EQ( s1.get(), 0.0 );

	ScalarWithUnits s2(100.0,UNITS_TEMPERATURE_CELSIUS);
	ASSERT_EQ( s2.get_units(), UNITS_TEMPERATURE_CELSIUS );
	ASSERT_DOUBLE_EQ( s2.get(), 100.0 );

	ScalarWithUnits s3(100.0,"degC");
	ASSERT_EQ( s3.get_units(), UNITS_TEMPERATURE_CELSIUS );
	ASSERT_DOUBLE_EQ( s3.get(), 100.0 );

	// copy constructor
	ScalarWithUnits s4( s1 );
	ASSERT_EQ( s4.get_units(), UNITS_NONE );
	ASSERT_DOUBLE_EQ( s4.get(), 0.0 );
	s1.set( 50.0, UNITS_TEMPERATURE_CELSIUS );
	ASSERT_EQ( s1.get_units(), UNITS_TEMPERATURE_CELSIUS );
	ASSERT_DOUBLE_EQ( s1.get(), 50.0 );
	ASSERT_EQ( s4.get_units(), UNITS_NONE );
	ASSERT_DOUBLE_EQ( s4.get(), 0.0 );

	// move constructor?
	ScalarWithUnits s5( std::move(ScalarWithUnits(100.0,UNITS_TEMPERATURE_CELSIUS)) );
	ASSERT_EQ( s5.get_units(), UNITS_TEMPERATURE_CELSIUS );
	ASSERT_DOUBLE_EQ( s5.get(), 100.0 );
}

TEST(ScalarWithUnitsTest, GetTests) {
	ScalarWithUnits freezing( 0.0, UNITS_TEMPERATURE_CELSIUS );
	ASSERT_DOUBLE_EQ( freezing.get(), 0.0 );
	ASSERT_DOUBLE_EQ( freezing.get_as(UNITS_TEMPERATURE_FAHRENHEIT), 32.0 );
	ASSERT_DOUBLE_EQ( freezing.get_as("Kelvin"), 273.15 );
	ASSERT_EQ( freezing.get_units(), UNITS_TEMPERATURE_CELSIUS);
}

TEST(ScalarWithUnitsTest, SetTests) {
	ScalarWithUnits freezing( 0.0, UNITS_TEMPERATURE_CELSIUS );

	// set temperature to 212 without changing units
	freezing.set_value( 212.0 );
	ASSERT_DOUBLE_EQ( freezing.get(), 212.0 );

	// set units to Fahrenheit without changing value
	freezing.set_units(UNITS_TEMPERATURE_FAHRENHEIT);
	ASSERT_DOUBLE_EQ( freezing.get_as(UNITS_TEMPERATURE_CELSIUS), 100.0 );

	// Set units to Kelvin using string
	freezing.set_units("Kelvin");
	ASSERT_DOUBLE_EQ( freezing.get_as(UNITS_TEMPERATURE_CELSIUS), -61.15 );

	// Set units to 32 degrees fahrenheit
	freezing.set( 32.0, UNITS_TEMPERATURE_FAHRENHEIT );
	ASSERT_DOUBLE_EQ( freezing.get_as(UNITS_TEMPERATURE_CELSIUS), 0.0 );

	// Set units to 100 degrees celsius
	freezing.set( 100.0, "c" );
	ASSERT_DOUBLE_EQ( freezing.get_as(UNITS_TEMPERATURE_FAHRENHEIT), 212.0 );
}

TEST(ScalarWithUnitsTest, ConvertTests) {
	ScalarWithUnits freezing( 0.0, UNITS_TEMPERATURE_CELSIUS );
	freezing.convert_units(UNITS_TEMPERATURE_FAHRENHEIT);
	ASSERT_DOUBLE_EQ(freezing.get(),32.0);
	freezing.convert_units("K");
	ASSERT_DOUBLE_EQ(freezing.get(),273.15);
}

TEST(ScalarWithUnitsTest, OperatorTests) {
	ScalarWithUnits firstLeg( 10.0, UNITS_DISTANCE_METERS );
	ScalarWithUnits secondLeg( 20.0, UNITS_DISTANCE_METERS );
	ScalarWithUnits thirdLeg( 0.03, UNITS_DISTANCE_KILOMETERS );

	ScalarWithUnits firstAndSecond = firstLeg + secondLeg;
	ScalarWithUnits secondAndThird = secondLeg + thirdLeg;
	ScalarWithUnits secondMinusThird = secondLeg - thirdLeg;

	ASSERT_DOUBLE_EQ( firstAndSecond.get(), 30.0 );
	ASSERT_DOUBLE_EQ( secondAndThird.get(), 50.0 );
	ASSERT_DOUBLE_EQ( secondMinusThird.get(), -10.0 );
	ASSERT_EQ( secondMinusThird.get_units(), UNITS_DISTANCE_METERS );

	firstLeg += secondLeg;
	thirdLeg -= secondLeg;
	ASSERT_DOUBLE_EQ( firstLeg.get(), 30.0 );
	ASSERT_DOUBLE_EQ( thirdLeg.get(), 0.01 );
}

