#include "NCPAUnits.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace NCPA;
using namespace testing;

TEST(unitsTest, EmptyStringConvertsToUnitsNone) {
	ASSERT_EQ( Units::fromString(""), units_t::NONE );
}

TEST(unitsTest, EqualityOperatorWorks) {
	ASSERT_TRUE( units_t::DISTANCE_METERS == units_t::DISTANCE_METERS );
	ASSERT_FALSE( units_t::DISTANCE_METERS == units_t::NONE );
}

TEST(unitsTest, InequalityOperatorWorks) {
	ASSERT_FALSE( units_t::DISTANCE_METERS != units_t::DISTANCE_METERS );
	ASSERT_TRUE( units_t::DISTANCE_METERS != units_t::NONE );
}

TEST(unitsTest,EqualOrNoneFunctionTrueForSameUnits) {
	ASSERT_TRUE( equal_or_none( units_t::DISTANCE_METERS, units_t::DISTANCE_METERS ) );
}

TEST(unitsTest,EqualOrNoneFunctionTrueForOneUnitIsNone) {
	ASSERT_TRUE( equal_or_none( units_t::NONE, units_t::DISTANCE_METERS ) );
	ASSERT_TRUE( equal_or_none( units_t::TEMPERATURE_KELVIN, units_t::NONE ) );
}

TEST(unitsTest,EqualOrNoneFunctionFalseForDifferentActualUnits) {
	ASSERT_FALSE( equal_or_none( units_t::DISTANCE_METERS, units_t::TEMPERATURE_KELVIN ) );
}

TEST(unitsTest, ConvertDistanceUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DISTANCE_METERS,units_t::DISTANCE_METERS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DISTANCE_KILOMETERS,units_t::DISTANCE_KILOMETERS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DISTANCE_METERS,units_t::DISTANCE_KILOMETERS), 0.001);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DISTANCE_KILOMETERS,units_t::DISTANCE_METERS), 1000.0);
}

TEST(unitsTest, ConvertTemperatureUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::TEMPERATURE_KELVIN,units_t::TEMPERATURE_KELVIN), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::TEMPERATURE_CELSIUS,units_t::TEMPERATURE_CELSIUS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::TEMPERATURE_FAHRENHEIT,units_t::TEMPERATURE_FAHRENHEIT), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::TEMPERATURE_KELVIN,units_t::TEMPERATURE_CELSIUS), -272.15);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::TEMPERATURE_CELSIUS,units_t::TEMPERATURE_KELVIN), 274.15);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::TEMPERATURE_KELVIN,units_t::TEMPERATURE_FAHRENHEIT), -457.87);
	ASSERT_DOUBLE_EQ( Units::convert(41.0,units_t::TEMPERATURE_FAHRENHEIT,units_t::TEMPERATURE_KELVIN), 278.15);
	ASSERT_DOUBLE_EQ( Units::convert(100.0,units_t::TEMPERATURE_CELSIUS,units_t::TEMPERATURE_FAHRENHEIT), 212.0);
	ASSERT_DOUBLE_EQ( Units::convert(41.0,units_t::TEMPERATURE_FAHRENHEIT,units_t::TEMPERATURE_CELSIUS), 5.0);
}

TEST(unitsTest, ConvertSpeedUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::SPEED_METERS_PER_SECOND,units_t::SPEED_METERS_PER_SECOND), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::SPEED_KILOMETERS_PER_SECOND,units_t::SPEED_KILOMETERS_PER_SECOND), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::SPEED_METERS_PER_SECOND,units_t::SPEED_KILOMETERS_PER_SECOND), 0.001);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::SPEED_KILOMETERS_PER_SECOND,units_t::SPEED_METERS_PER_SECOND), 1000.0);
}

TEST(unitsTest, ConvertPressureUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_PASCALS,units_t::PRESSURE_PASCALS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_HECTOPASCALS,units_t::PRESSURE_HECTOPASCALS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_MILLIBARS,units_t::PRESSURE_MILLIBARS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_ATMOSPHERES,units_t::PRESSURE_ATMOSPHERES), 1.0);

	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_PASCALS,units_t::PRESSURE_HECTOPASCALS), 0.01);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_PASCALS,units_t::PRESSURE_MILLIBARS), 0.01);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_PASCALS,units_t::PRESSURE_ATMOSPHERES), 0.00000986923);

	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_HECTOPASCALS,units_t::PRESSURE_PASCALS), 100.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_HECTOPASCALS,units_t::PRESSURE_MILLIBARS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_HECTOPASCALS,units_t::PRESSURE_ATMOSPHERES), 0.000986923);

	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_MILLIBARS,units_t::PRESSURE_HECTOPASCALS), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_MILLIBARS,units_t::PRESSURE_PASCALS), 100.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_MILLIBARS,units_t::PRESSURE_ATMOSPHERES), 0.000986923);

	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_ATMOSPHERES,units_t::PRESSURE_HECTOPASCALS), 1013.25);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_ATMOSPHERES,units_t::PRESSURE_PASCALS), 101325.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::PRESSURE_ATMOSPHERES,units_t::PRESSURE_MILLIBARS), 1013.25);
}

TEST(unitsTest, ConvertDensityUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER,units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER,units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER,units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER), 0.001);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER,units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER), 1000.0);
}

TEST(unitsTest, ConvertDirectionUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH,units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(1.0,units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST,units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST), 1.0);
	ASSERT_DOUBLE_EQ( Units::convert(0.0,units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH,units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST), 90.0);
	ASSERT_DOUBLE_EQ( Units::convert(180.0,units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST,units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH), 270.0);
}

TEST(unitsTest, ConvertAngleUnits) {
	ASSERT_DOUBLE_EQ( Units::convert(180.0,units_t::ANGLE_DEGREES,units_t::ANGLE_DEGREES), 180.0);
	ASSERT_DOUBLE_EQ( Units::convert(2.0,units_t::ANGLE_RADIANS,units_t::ANGLE_RADIANS), 2.0);
	ASSERT_DOUBLE_EQ( Units::convert(180.0,units_t::ANGLE_DEGREES,units_t::ANGLE_RADIANS), M_PI);
	ASSERT_DOUBLE_EQ( Units::convert(3.0*M_PI/2.0,units_t::ANGLE_RADIANS,units_t::ANGLE_DEGREES), 270.0);
}

TEST(unitsTest, ConvertAttenuationUnits) {
	ASSERT_DOUBLE_EQ( Units::convert( 10.0, units_t::ATTENUATION_NEPERS_PER_METER, units_t::ATTENUATION_DECIBELS_PER_KILOMETER ),  86858.89638);
	ASSERT_DOUBLE_EQ( Units::convert( 10.0, units_t::ATTENUATION_NEPERS_PER_METER, units_t::ATTENUATION_DECIBELS_PER_METER ),  86.85889638);
	ASSERT_DOUBLE_EQ( Units::convert( 10.0, units_t::ATTENUATION_DECIBELS_PER_KILOMETER, units_t::ATTENUATION_NEPERS_PER_METER ),  0.00115129255 );
	ASSERT_DOUBLE_EQ( Units::convert( 10.0, units_t::ATTENUATION_DECIBELS_PER_KILOMETER, units_t::ATTENUATION_DECIBELS_PER_METER ),  0.01 );
	ASSERT_DOUBLE_EQ( Units::convert( 10.0, units_t::ATTENUATION_DECIBELS_PER_METER, units_t::ATTENUATION_NEPERS_PER_METER ),  1.15129255 );
	ASSERT_DOUBLE_EQ( Units::convert( 10.0, units_t::ATTENUATION_DECIBELS_PER_METER, units_t::ATTENUATION_DECIBELS_PER_KILOMETER ),  10000.0 );
}

TEST(unitsTest, ConvertRatioUnits) {
	ASSERT_DOUBLE_EQ( Units::convert( 0.15, units_t::RATIO_DECIMAL, units_t::RATIO_PERCENT ),  15);
	ASSERT_DOUBLE_EQ( Units::convert( 0.15, units_t::RATIO_PERCENT, units_t::RATIO_DECIMAL ),  0.0015);
}

TEST(unitsTest, ConvertArrayTest) {
	double in[5] = { 1.0, 2.0, 4.0, -2.0, -5.0 };
	double out[5];
	Units::convert( in, 5, "km", "m", out );
	ASSERT_DOUBLE_EQ( out[0], 1000.0 );
	ASSERT_DOUBLE_EQ( out[1], 2000.0 );
	ASSERT_DOUBLE_EQ( out[2], 4000.0 );
	ASSERT_DOUBLE_EQ( out[3], -2000.0 );
	ASSERT_DOUBLE_EQ( out[4], -5000.0 );
}

TEST(unitsTest, ConvertUsingStringsTest) {
	double km = 5.0;
	ASSERT_DOUBLE_EQ( Units::convert( km, "km", "m" ), 5000.0 );
	ASSERT_DOUBLE_EQ( Units::convert( km, units_t::DISTANCE_KILOMETERS, "m" ), 5000.0 );
	ASSERT_DOUBLE_EQ( Units::convert( km, "km", units_t::DISTANCE_METERS ), 5000.0 );
	double m;
	Units::convert( &km, 1, "km", "m", &m );
	ASSERT_DOUBLE_EQ( m, 5000.0 );
	m=0;
	Units::convert( &km, 1, units_t::DISTANCE_METERS, "km", &m );
	ASSERT_DOUBLE_EQ( m, 0.005 );
	m=0;
	Units::convert( &km, 1, "km", units_t::DISTANCE_METERS, &m );
	ASSERT_DOUBLE_EQ( m, 5000.0 );

}

TEST(unitsTest, FromStringTest) {
	ASSERT_EQ( Units::fromString( "" ), units_t::NONE );
	ASSERT_EQ( Units::fromString( "N/A" ), units_t::NONE );
	ASSERT_EQ( Units::fromString( "K" ), units_t::TEMPERATURE_KELVIN );
	ASSERT_EQ( Units::fromString( "DEGK" ), units_t::TEMPERATURE_KELVIN );
	ASSERT_EQ( Units::fromString( "DEG K" ), units_t::TEMPERATURE_KELVIN );
	ASSERT_EQ( Units::fromString( "KELVIN" ), units_t::TEMPERATURE_KELVIN );
	ASSERT_EQ( Units::fromString( "C" ), units_t::TEMPERATURE_CELSIUS );
	ASSERT_EQ( Units::fromString( "DEGC" ), units_t::TEMPERATURE_CELSIUS );
	ASSERT_EQ( Units::fromString( "DEG C" ), units_t::TEMPERATURE_CELSIUS );
	ASSERT_EQ( Units::fromString( "DEGREES C" ), units_t::TEMPERATURE_CELSIUS );
	ASSERT_EQ( Units::fromString( "CELSIUS" ), units_t::TEMPERATURE_CELSIUS );
	ASSERT_EQ( Units::fromString( "CENTIGRADE" ), units_t::TEMPERATURE_CELSIUS );
	ASSERT_EQ( Units::fromString( "F" ), units_t::TEMPERATURE_FAHRENHEIT );
	ASSERT_EQ( Units::fromString( "DEGF" ), units_t::TEMPERATURE_FAHRENHEIT );
	ASSERT_EQ( Units::fromString( "DEG F" ), units_t::TEMPERATURE_FAHRENHEIT );
	ASSERT_EQ( Units::fromString( "DEGREES F" ), units_t::TEMPERATURE_FAHRENHEIT );
	ASSERT_EQ( Units::fromString( "FAHRENHEIT" ), units_t::TEMPERATURE_FAHRENHEIT );
	ASSERT_EQ( Units::fromString( "KM" ), units_t::DISTANCE_KILOMETERS );
	ASSERT_EQ( Units::fromString( "kilometers" ), units_t::DISTANCE_KILOMETERS );
	ASSERT_EQ( Units::fromString( "m" ), units_t::DISTANCE_METERS );
	ASSERT_EQ( Units::fromString( "meters" ), units_t::DISTANCE_METERS );
	ASSERT_EQ( Units::fromString( "m/s" ), units_t::SPEED_METERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "mps" ), units_t::SPEED_METERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "mpers" ), units_t::SPEED_METERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "m per s" ), units_t::SPEED_METERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "meters per second" ), units_t::SPEED_METERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "km/s" ), units_t::SPEED_KILOMETERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "kmps" ), units_t::SPEED_KILOMETERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "kmpers" ), units_t::SPEED_KILOMETERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "km per s" ), units_t::SPEED_KILOMETERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "kilometers per second" ), units_t::SPEED_KILOMETERS_PER_SECOND );
	ASSERT_EQ( Units::fromString( "Pa" ), units_t::PRESSURE_PASCALS );
	ASSERT_EQ( Units::fromString( "Pascal" ), units_t::PRESSURE_PASCALS );
	ASSERT_EQ( Units::fromString( "Pascals" ), units_t::PRESSURE_PASCALS );
	ASSERT_EQ( Units::fromString( "mbar" ), units_t::PRESSURE_MILLIBARS );
	ASSERT_EQ( Units::fromString( "millibar" ), units_t::PRESSURE_MILLIBARS );
	ASSERT_EQ( Units::fromString( "millibars" ), units_t::PRESSURE_MILLIBARS );
	ASSERT_EQ( Units::fromString( "hPa" ), units_t::PRESSURE_HECTOPASCALS );
	ASSERT_EQ( Units::fromString( "hectopascals" ), units_t::PRESSURE_HECTOPASCALS );
	ASSERT_EQ( Units::fromString( "atm" ), units_t::PRESSURE_ATMOSPHERES );
	ASSERT_EQ( Units::fromString( "atmospheres" ), units_t::PRESSURE_ATMOSPHERES );
	ASSERT_EQ( Units::fromString( "kg/m3" ), units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER );
	ASSERT_EQ( Units::fromString( "kgpm3" ), units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER );
	ASSERT_EQ( Units::fromString( "KILOGRAMS PER CUBIC METER" ), units_t::DENSITY_KILOGRAMS_PER_CUBIC_METER );
	ASSERT_EQ( Units::fromString( "g/cm3" ), units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER );
	ASSERT_EQ( Units::fromString( "gpcm3" ), units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER );
	ASSERT_EQ( Units::fromString( "GRAMS PER CUBIC CENTIMETER" ), units_t::DENSITY_GRAMS_PER_CUBIC_CENTIMETER );
	ASSERT_EQ( Units::fromString( "degrees clockwise from north" ), units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
	ASSERT_EQ( Units::fromString( "deg cw from n" ), units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
	ASSERT_EQ( Units::fromString( "azimuth" ), units_t::DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
	ASSERT_EQ( Units::fromString( "degrees counterclockwise from east" ), units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST );
	ASSERT_EQ( Units::fromString( "deg ccw from e" ), units_t::DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST );
	ASSERT_EQ( Units::fromString( "deg" ), units_t::ANGLE_DEGREES );
	ASSERT_EQ( Units::fromString( "degrees" ), units_t::ANGLE_DEGREES );
	ASSERT_EQ( Units::fromString( "rad" ), units_t::ANGLE_RADIANS );
	ASSERT_EQ( Units::fromString( "radians" ), units_t::ANGLE_RADIANS );
}

TEST(unitsTest, BadConversionTest) {
	double in = 1.0, out;
	EXPECT_THROW( {out = Units::convert( in, units_t::PRESSURE_MILLIBARS, units_t::ANGLE_DEGREES );},
			NCPA::invalid_conversion );
}

TEST(unitsTest, CanConvertReturnsTrueBetweenCompatibleUnits) {
	ASSERT_TRUE( Units::can_convert(
			units_t::DISTANCE_METERS, units_t::DISTANCE_KILOMETERS ) );
	ASSERT_TRUE( Units::can_convert( "m", units_t::DISTANCE_KILOMETERS ) );
	ASSERT_TRUE( Units::can_convert( units_t::DISTANCE_METERS, "km" ) );
	ASSERT_TRUE( Units::can_convert( "m", "km" ) );
}

TEST(unitsTest, CanConvertReturnsFalseBetweenIncompatibleUnits) {
	ASSERT_FALSE( Units::can_convert(
			units_t::DISTANCE_METERS, units_t::TEMPERATURE_KELVIN ) );
	ASSERT_FALSE( Units::can_convert( "m", units_t::TEMPERATURE_KELVIN ) );
	ASSERT_FALSE( Units::can_convert( units_t::DISTANCE_METERS, "K" ) );
	ASSERT_FALSE( Units::can_convert( "m", "K" ) );
}



