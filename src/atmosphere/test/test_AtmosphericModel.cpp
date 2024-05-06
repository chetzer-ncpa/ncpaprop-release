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



TEST(AtmosphericModelTest,SoundspeedFromTemperatureIsCorrect) {
	double T = 273.15;  // 0 celsius
	units_t T_units = Units::fromString("K");
	units_t C_units = Units::fromString("m/s");
	ASSERT_NEAR( AtmosphericModel::soundspeed_from_temperature(T, T_units, C_units), 331.28, 1.0e-2 );
	ASSERT_NEAR( AtmosphericModel::soundspeed_from_temperature(0.0, Units::fromString("C"), C_units), 331.287895946, 1.0e-6 );
	T = 323.15;
	ASSERT_NEAR( AtmosphericModel::soundspeed_from_temperature(T, T_units, C_units), 360.34, 1.0e-2 );
	ASSERT_NEAR( AtmosphericModel::soundspeed_from_temperature(T, "K", "m/s"), 360.34, 1.0e-2 );
}

TEST(AtmosphericModelTest,TemperatureFromSoundspeedIsCorrect) {
	double C = 331.2879;  // 0 celsius
	units_t T_units = Units::fromString("C");
	units_t C_units = Units::fromString("m/s");
	ASSERT_NEAR( AtmosphericModel::temperature_from_soundspeed(C, C_units, T_units), 0.0, 1.0e-2 );
	ASSERT_NEAR( AtmosphericModel::temperature_from_soundspeed(C, "m/s", "C"), 0.0, 1.0e-2 );

}

TEST(AtmosphericModelTest,SoundspeedFromPressureAndDensityIsCorrect) {
	double P = 1.0;
	units_t P_units = Units::fromString("atmospheres");
	double D = 1.225;
	units_t D_units = Units::fromString("kg/m3");
	units_t C_units = Units::fromString("m/s");
	ASSERT_NEAR( AtmosphericModel::soundspeed_from_pressure_density(
					P, P_units, D, D_units, C_units
				), 340.293, 1.0e-2 );
	ASSERT_NEAR( AtmosphericModel::soundspeed_from_pressure_density(
					P, "atm", D, "kg/m3", "m/s"
				), 340.293, 1.0e-2 );
}

// 1.225 -> 0.001225
TEST(AtmosphericModelTest,DensityFromTemperatureAndPressureIsCorrect) {
	double P = 1.0;
	units_t P_units = Units::fromString("atmospheres");
	double T = 15.0;
	units_t T_units = Units::fromString("C");
	units_t D_units = Units::fromString("g/cm3");
	ASSERT_NEAR( AtmosphericModel::density_from_temperature_pressure(
					T, T_units, P, P_units, D_units
			), 0.001225, 1e-5);
	ASSERT_NEAR( AtmosphericModel::density_from_temperature_pressure(
					T, "C", P, "atm", "g/cm3"
			), 0.001225, 1e-5);
}
