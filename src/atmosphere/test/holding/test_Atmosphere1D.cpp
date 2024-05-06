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

// test fixture
class Atmosphere1DTest : public ::testing::Test {
 protected:

	// Initializations and other setup
	void SetUp() override {
		// test constructors
		a1 = Atmosphere1D( 10, alts, Units::fromString("km") );
		a2 = Atmosphere1D( 10, alts, "km" );

		// test copy constructor
		altv = VectorWithUnits(10,alts,"km");
		uv = VectorWithUnits(10,uwinds,"m/s");
		vv = VectorWithUnits(10,vwinds,"m/s");
		tv = VectorWithUnits(10,temps,"C");
		cv = VectorWithUnits(10,speeds,"m/s");
		pv = VectorWithUnits(10,pressures,"atmospheres");
		rv = VectorWithUnits(10,densities,"kg/m3");
		z0 = ScalarWithUnits(54.0,"m");
		a3 = a1;
		a3.add_property( "U", uv );
		a3.add_property( "T", tv );
		a3.add_property( "V", vv );
		a3.add_property( "C0", cv );
		a3.add_property( "P", pv );
		a3.add_property( "RHO", rv );
		a3.add_property( "Z0", z0 );

		// read profiles from files
		a4 = Atmosphere1D(testfilename_profile);
		a5 = Atmosphere1D(testfilename_profile_noheader, testfilename_header);
	}

	// If there's any cleanup other than normal destructors
	//	void TearDown() override {}

	// class members here
	Atmosphere1D a0, a1, a2, a3, a4, a5;
	double temps[10] = { 0.0, 12.0, 8.0, 3.0, 0.0, -1.0, -1.3, -2.0, -3.0, -5.0 };
	double uwinds[10] = { 0, 4, 5, 2, -1, -7, -2, 2, 5, 10 };
	double vwinds[10] = { 1, 4, 5, 2, -1, -7, -2, 2, 5, 10 };
	double alts[10]  = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double speeds[10] = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	double pressures[10] = { 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 };
	double densities[10] = {1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5};
	VectorWithUnits altv, tv, uv, vv, cv, pv, rv;
	ScalarWithUnits z0;

	string  testfilename_profile = "testprofile.dat",
			testfilename_profile_noheader = "testprofile_noheader.dat",
			testfilename_header = "testheader.dat";
};

TEST_F(Atmosphere1DTest,DirectConstructorsWork) {
	ASSERT_EQ(a0.nz(),0);
	ASSERT_EQ(a1.nz(),10);
	ASSERT_EQ(a2.nz(),10);
	ASSERT_EQ(a3.nz(),10);
}

TEST_F(Atmosphere1DTest,FileConstructorsWork) {
	ASSERT_TRUE(a4.contains_scalar("Z0"));
	ASSERT_TRUE(a4.contains_vector("U"));
	ASSERT_TRUE(a4.contains_vector("V"));
	ASSERT_TRUE(a4.contains_vector("W"));
	ASSERT_TRUE(a4.contains_vector("T"));
	ASSERT_TRUE(a4.contains_vector("RHO"));
	ASSERT_TRUE(a4.contains_vector("P"));

	ASSERT_TRUE(a5.contains_scalar("Z0"));
	ASSERT_TRUE(a5.contains_vector("U"));
	ASSERT_TRUE(a5.contains_vector("V"));
	ASSERT_TRUE(a5.contains_vector("W"));
	ASSERT_TRUE(a5.contains_vector("T"));
	ASSERT_TRUE(a5.contains_vector("RHO"));
	ASSERT_TRUE(a5.contains_vector("P"));

	ASSERT_DOUBLE_EQ( a4.get("Z0"), 0.0 );
	ASSERT_DOUBLE_EQ( a5.get("Z0"), 0.0 );

	// check a couple of random points for each parameter
	double ztest = 15.5;
	ASSERT_DOUBLE_EQ(a4.get("U",ztest), 0.0933147);
	ASSERT_DOUBLE_EQ(a4.get("V",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a4.get("W",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a4.get("T",ztest), 223.9038);
	ASSERT_DOUBLE_EQ(a4.get("RHO",ztest), 1.76811e-4);
	ASSERT_DOUBLE_EQ(a4.get("P",ztest), 113.6392);

	ASSERT_DOUBLE_EQ(a5.get("U",ztest), 0.0933147);
	ASSERT_DOUBLE_EQ(a5.get("V",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a5.get("W",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a5.get("T",ztest), 223.9038);
	ASSERT_DOUBLE_EQ(a5.get("RHO",ztest), 1.76811e-4);
	ASSERT_DOUBLE_EQ(a5.get("P",ztest), 113.6392);

	ztest = 85.9;
	ASSERT_DOUBLE_EQ(a4.get("U",ztest), 6.71232);
	ASSERT_DOUBLE_EQ(a4.get("V",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a4.get("W",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a4.get("T",ztest), 128.4423);
	ASSERT_DOUBLE_EQ(a4.get("RHO",ztest), 7.77123e-9);
	ASSERT_DOUBLE_EQ(a4.get("P",ztest), 2.865203e-3);

	ASSERT_DOUBLE_EQ(a5.get("U",ztest), 6.71232);
	ASSERT_DOUBLE_EQ(a5.get("V",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a5.get("W",ztest), 0.0);
	ASSERT_DOUBLE_EQ(a5.get("T",ztest), 128.4423);
	ASSERT_DOUBLE_EQ(a5.get("RHO",ztest), 7.77123e-9);
	ASSERT_DOUBLE_EQ(a5.get("P",ztest), 2.865203e-3);

}

TEST_F(Atmosphere1DTest,AddPropertyWorksWithVectorWithUnits) {
	a1.add_property("T",tv);
	ASSERT_DOUBLE_EQ(a1.get("T",1.0),12.0);
}

TEST_F(Atmosphere1DTest,AddPropertyWorksWithArrayAndUnits) {
	a1.add_property("U",10,uwinds,Units::fromString("m/s"));
	ASSERT_DOUBLE_EQ(a1.get("U",5.0),-7.0);
}

TEST_F(Atmosphere1DTest,AddPropertyWorksWithArrayAndString) {
	a1.add_property("U",10,uwinds,"m/s");
	ASSERT_DOUBLE_EQ(a1.get("U",5.0),-7.0);
}

TEST_F(Atmosphere1DTest,AddPropertyWorksWithScalarWithUnits) {
	a1.add_property( "Z0", ScalarWithUnits(54.0,"m") );
	ASSERT_TRUE(a1.contains_key("Z0"));
	ASSERT_DOUBLE_EQ(a1.get("Z0"),54.0);
}

TEST_F(Atmosphere1DTest,AddPropertyWorksWithDoubleAndUnits) {
	a1.add_property( "Z0", 54.0, Units::fromString("m") );
	ASSERT_TRUE(a1.contains_key("Z0"));
	ASSERT_DOUBLE_EQ(a1.get("Z0"),54.0);
}

TEST_F(Atmosphere1DTest,AddPropertyWorksWithDoubleAndString) {
	a1.add_property( "Z0", 54.0, "m" );
	ASSERT_TRUE(a1.contains_key("Z0"));
	ASSERT_DOUBLE_EQ(a1.get("Z0"),54.0);
}

TEST_F(Atmosphere1DTest,CopyConstructorWorks) {
	a1.add_property("U",10,uwinds,"m/s");
	Atmosphere1D a5( a1 );
	ASSERT_EQ(a5.nz(),10);
	ASSERT_DOUBLE_EQ(a5.get("U",5.0),-7.0);
}

TEST_F(Atmosphere1DTest,SwapWorks) {
	swap(a0,a1);
	ASSERT_EQ(a0.nz(),10);
	ASSERT_EQ(a1.nz(),0);
}

TEST_F(Atmosphere1DTest,AssignmentOperatorWorks) {
	a0 = a1;
	ASSERT_EQ(a0.nz(),10);
	ASSERT_EQ(a1.nz(),10);
}

TEST_F(Atmosphere1DTest,ResetAltitudeVectorWorks) {
	a1.add_property("U",10,uwinds,"m/s");
	double newz[5] = { 0, 2, 4, 6, 8 };
	VectorWithUnits v( 5, newz, "km" );
	a1.reset_altitude_vector(v);
	ASSERT_EQ(a1.nz(),5);
	ASSERT_FALSE( a1.contains_key("U") );
}

TEST_F(Atmosphere1DTest,ResetAltitudeVectorWorksWithArrayAndUnits) {
	a1.add_property("U",10,uwinds,"m/s");
	double newz[5] = { 0, 2, 4, 6, 8 };
	a1.reset_altitude_vector(5, newz, Units::fromString("km") );
	ASSERT_EQ(a1.nz(),5);
	ASSERT_FALSE( a1.contains_key("U") );
}

TEST_F(Atmosphere1DTest,ResetAltitudeVectorWorksWithArrayAndString) {
	a1.add_property("U",10,uwinds,"m/s");
	double newz[5] = { 0, 2, 4, 6, 8 };
	a1.reset_altitude_vector(5, newz, "km" );
	ASSERT_EQ(a1.nz(),5);
	ASSERT_FALSE( a1.contains_key("U") );
}

TEST_F(Atmosphere1DTest,CopyVectorPropertyWorks) {
	a3.copy_vector_property( "U", "UCOPY" );
	ASSERT_DOUBLE_EQ(a3.get("U",3.4),a3.get("UCOPY",3.4));
	ASSERT_DOUBLE_EQ(a3.get("U",8.1),a3.get("UCOPY",8.1));
	ASSERT_EQ(a3.get_property_units("U"),a3.get_property_units("UCOPY"));
}

TEST_F(Atmosphere1DTest,CopyScalarPropertyWorks) {
	a3.copy_scalar_property( "Z0", "Z0COPY" );
	ASSERT_DOUBLE_EQ(a3.get("Z0"),a3.get("Z0COPY"));
	ASSERT_EQ(a3.get_property_units("Z0"),a3.get_property_units("Z0COPY"));
}

TEST_F(Atmosphere1DTest,GetVectorPropertyObjectWorks) {
	NCPA::AtmosphericProperty1D *p = a3.get_vector_property_object("U");
	for (double d = 0.0; d <= 9.0; d += 0.1) {
		ASSERT_DOUBLE_EQ(p->get(d),a3.get("U",d));
	}
}

TEST_F(Atmosphere1DTest,GetConstVectorPropertyObjectWorks) {
	const NCPA::AtmosphericProperty1D *p = a3.get_const_vector_property_object("U");
	for (double d = 0.0; d <= 9.0; d += 0.1) {
		ASSERT_DOUBLE_EQ(p->get(d),a3.get("U",d));
	}
}

TEST_F(Atmosphere1DTest,GetScalarPropertyObjectWorks) {
	ScalarWithUnits *u = a3.get_scalar_property_object("Z0");
	ASSERT_DOUBLE_EQ(u->get(),54.0);
	u->set(25.0,"m");
	ASSERT_DOUBLE_EQ(a3.get("Z0"),25.0);
}

TEST_F(Atmosphere1DTest,GetConstScalarPropertyObjectWorks) {
	const ScalarWithUnits *u = a3.get_const_scalar_property_object("Z0");
	ASSERT_DOUBLE_EQ(u->get(),54.0);
}

TEST_F(Atmosphere1DTest,ResetSplinesDoesntAffectValues) {
	vector<double> a, u;
	for (double d = 0.0; d <= 9.0; d += 0.1) {
		a.push_back(d);
		u.push_back( a3.get("U",d) );
	}
	a3.reset_splines();
	auto ai = a.begin();
	auto ui = u.begin();
	for ( ; ai != a.end() && ui != u.end(); ++ai, ++ui) {
		ASSERT_DOUBLE_EQ( a3.get("U",*ai), *ui );
	}
}

TEST_F(Atmosphere1DTest,ResetSplinesAfterChangeWorksProperly) {
	AtmosphericProperty1D *p = a3.get_vector_property_object("U");
	p->second[5] = -10.0;
	a3.reset_splines();
	ASSERT_DOUBLE_EQ(a3.get("U",5.0),-10.0);
}

TEST_F(Atmosphere1DTest,RemovePropertyWorks) {
	a3.remove_property("U");
	ASSERT_FALSE(a3.contains_key("U"));
}

TEST_F(Atmosphere1DTest,GetAltitudeVectorWorks) {
	double buffer[10];
	a3.get_altitude_vector( buffer );
	for (size_t i = 0; i < 10; i++) {
		ASSERT_DOUBLE_EQ(buffer[i],alts[i]);
	}
}

TEST_F(Atmosphere1DTest,GetAltitudeVectorWithUnitsWorks) {
	double buffer[10];
	units_t u;
	a3.get_altitude_vector( buffer, &u );
	for (size_t i = 0; i < 10; i++) {
		ASSERT_DOUBLE_EQ(buffer[i],alts[i]);
	}
	ASSERT_EQ(u,Units::fromString("km"));
}

TEST_F(Atmosphere1DTest,GetPropertyVectorWorks) {
	double buffer[10];
	a3.get_property_vector( "U", buffer );
	for (size_t i = 0; i < 10; i++) {
		ASSERT_DOUBLE_EQ(buffer[i],uwinds[i]);
	}
}

TEST_F(Atmosphere1DTest,GetPropertyVectorWithUnitsWorks) {
	double buffer[10];
	units_t u;
	a3.get_property_vector( "U", buffer, &u );
	for (size_t i = 0; i < 10; i++) {
		ASSERT_DOUBLE_EQ(buffer[i],uwinds[i]);
	}
	ASSERT_EQ(u,Units::fromString("m/s"));
}

TEST_F(Atmosphere1DTest,GetThrowsOutOfRangeWithBadKey) {
	EXPECT_THROW( {double d = a3.get("W",3.0);},out_of_range);
}

TEST_F(Atmosphere1DTest,AddPropertyThrowsRuntimeErrorWithExistingKey) {
	EXPECT_THROW( {a3.add_property("U",uv);},invalid_argument);
}

TEST_F(Atmosphere1DTest,GetMinimumAltitudeWorks) {
	EXPECT_DOUBLE_EQ(a2.get_minimum_altitude(),0.0);
}

TEST_F(Atmosphere1DTest,GetMaximumAltitudeWorks) {
	EXPECT_DOUBLE_EQ(a2.get_maximum_altitude(),9.0);
}

TEST_F(Atmosphere1DTest,SoundSpeedFromTemperatureWorks) {
	a3.calculate_sound_speed_from_temperature( "C", "T", "m/s" );
	double expected[10] = {331.288,338.487,336.104,333.102,331.288,
			330.681,330.499,330.073,329.464,328.242};
	double buffer[10];
	units_t u;
	a3.get_property_vector( "C", buffer, &u );
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer,expected,0.1);
	EXPECT_EQ(u,Units::fromString("m/s"));
}

TEST_F(Atmosphere1DTest,TemperatureFromSoundSpeedWorks) {
	a3.calculate_temperature_from_sound_speed( "T0", "C0", "K" );
	double expected[10];
	for (size_t i = 0; i < 10; i++) {
		expected[i] = (speeds[i]*speeds[i] / 287.0 / 1.4);
	}
	double buffer[10];
	units_t u;
	a3.get_property_vector( "T0", buffer, &u );
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer,expected,0.1);
	EXPECT_EQ(u,Units::fromString("K"));
}

TEST_F(Atmosphere1DTest,SoundSpeedFromPressureAndDensityWorks) {
	a3.calculate_sound_speed_from_pressure_and_density( "C", "P", "RHO", "m/s" );
	double expected[10];
	for (size_t i = 0; i < 10; i++) {
		double p_Pa = Units::convert(pressures[i],"atmospheres","Pa");
		double d_kg = densities[i];
		expected[i] = std::sqrt( 1.4 * p_Pa / d_kg );
	}
	double buffer[10];
	units_t u;
	a3.get_property_vector( "C", buffer, &u );
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer,expected,0.1);
	EXPECT_EQ(u,Units::fromString("m/s"));
}

TEST_F(Atmosphere1DTest,DensityFromTemperatureAndPressureWorks) {
	a3.calculate_density_from_temperature_and_pressure( "D0", "T", "P", "kg/m3" );
	double expected[10];
	for (size_t i = 0; i < 10; i++) {
		double p_Pa = Units::convert(pressures[i],"atmospheres","Pa");
		double t_k = Units::convert(temps[i],"C","K");
		expected[i] = p_Pa / (t_k * 287.0);
	}
	double buffer[10];
	units_t u;
	a3.get_property_vector( "D0", buffer, &u );
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer,expected,0.1);
	EXPECT_EQ(u,Units::fromString("kg/m3"));
}

TEST_F(Atmosphere1DTest,CalculateWindSpeedWorks) {
	a3.calculate_wind_speed( "WS", "U", "V" );
	double buffer[10];
	units_t u;
	a3.get_property_vector( "WS", buffer, &u );
	double expected[10];
	for (size_t i = 0; i < 10; i++) {
		expected[i] = sqrt(uwinds[i]*uwinds[i] + vwinds[i]*vwinds[i]);
	}
	EXPECT_DOUBLE_ARRAY_EQ(10,buffer,expected);
}

TEST_F(Atmosphere1DTest,CalculateWindDirectionWorks) {
	a3.calculate_wind_direction( "WD", "U", "V" );
	double buffer[10];
	units_t u;
	a3.get_property_vector( "WD", buffer, &u );
	double expected[10];
	for (size_t i = 0; i < 10; i++) {
		expected[i] = normalizeAzimuth(90.0 - Units::convert(atan2(vwinds[i],uwinds[i]),"rad","deg"));
	}
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer,expected,0.1);
}

TEST_F(Atmosphere1DTest,CalculateWindComponentWorks) {
	a3.calculate_wind_direction( "WD", "U", "V" );
	a3.calculate_wind_speed( "WS", "U", "V" );
	a3.calculate_wind_component( "W90", "WS", "WD", 90.0 );
	a3.calculate_wind_component( "W180", "WS", "WD", 180.0 );
	double buffer90[10], buffer180[10];
	units_t u;
	a3.get_property_vector( "W90", buffer90, &u );
	a3.get_property_vector( "W180", buffer180, &u );
//	a3.get_property_vector( "WS", speedbuffer, &u );
	double expected[10];
	for (size_t i = 0; i < 10; i++) {
		expected[i] = -vwinds[i];
	}
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer90,uwinds,0.1);
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer180,expected,0.1);
}

TEST_F(Atmosphere1DTest,CalculateEffectiveSoundSpeedWorks) {
	a3.calculate_wind_direction( "WD", "U", "V" );
	a3.calculate_wind_speed( "WS", "U", "V" );
	a3.calculate_wind_component( "W90", "WS", "WD", 90.0 );
	a3.calculate_wind_component( "W180", "WS", "WD", 180.0 );
	a3.calculate_sound_speed_from_temperature( "C", "T", "m/s" );
	a3.calculate_effective_sound_speed( "CEFF90", "C", "W90" );
	a3.calculate_effective_sound_speed( "CEFF180", "C", "W180" );
	double buffer90[10], buffer180[10];
	a3.get_property_vector("CEFF90", buffer90);
	a3.get_property_vector("CEFF180", buffer180);
	double expected90[10], expected180[10];
	for (size_t i = 0; i < 10; i++) {
		double cexp = sqrt( 287.0 * 1.4 * (273.15+temps[i]));
		expected90[i] = cexp + uwinds[i];
		expected180[i] = cexp - vwinds[i];
	}
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer90,expected90,0.001);
	EXPECT_DOUBLE_ARRAY_NEAR(10,buffer180,expected180,0.001);
}


// Still need to test:
// read_values_from_stream
// read_headers_from_stream
// print_atmosphere
