#include "ncpaprop_common.h"
#include "ncpaprop_atmosphere.h"


#include <vector>


// set up the atmosphere for normal NCPAprop usage (units, building
// derived properties, etc)
void NCPA::setup_ncpaprop_atmosphere( NCPA::Atmosphere1D *atm_profile ) {

	// altitude units in meters
	atm_profile->convert_altitude_units( Units::fromString( "m" ) );

	// if the user supplies the static sound speed, keep it.  Use it
	// to calculate the temperature if it's absent
	if (atm_profile->contains_vector( "C0" ) ) {
		// std::cout << "Using user-supplied static sound speed" << std::endl;
		atm_profile->copy_vector_property( "C0", "_C0_" );
		atm_profile->convert_property_units( "_C0_", Units::fromString( "m/s" ) );

		// do we need to calculate temperature?
		if (!atm_profile->contains_vector("T" )) {
			// std::cout << "No temperature provided, calculating from sound speed"
			// 		  << std::endl;
			atm_profile->calculate_temperature_from_sound_speed( "T", "_C0_",
				Units::fromString( "K" ) );
		}
	} else {
		if (atm_profile->contains_vector( "P" ) && atm_profile->contains_vector("RHO")) {
			// std::cout << "Calculating sound speed from pressure and density." << std::endl;
			atm_profile->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO",
				Units::fromString( "m/s" ) );

			if (!atm_profile->contains_vector("T" )) {
				// std::cout << "No temperature provided, calculating from sound speed"
				// 		  << std::endl;
				atm_profile->calculate_temperature_from_sound_speed( "T", "_C0_",
					Units::fromString( "K" ) );
			}
		} else if (atm_profile->contains_vector( "T" )) {
			// std::cout << "Calculating sound speed from temperature" << std::endl;
			atm_profile->calculate_sound_speed_from_temperature( "_C0_", "T",
				Units::fromString( "m/s" ) );
		} else {
			throw std::runtime_error( "Can't calculate static sound speed, either temperature T or pressure P and density RHO missing" );
		}
	}

	// set units
	if (atm_profile->contains_vector("U")) {
		atm_profile->convert_property_units( "U", Units::fromString( "m/s" ) );
	} else {
		throw std::runtime_error( "U wind vector not found in atmosphere file" );
	}
	if (atm_profile->contains_vector("V")) {
		atm_profile->convert_property_units( "V", Units::fromString( "m/s" ) );
	} else {
		throw std::runtime_error( "V wind vector not found in atmosphere file" );
	}
	if (atm_profile->contains_vector("T")) {
		atm_profile->convert_property_units( "T", Units::fromString( "K" ) );
	} else {
		throw std::runtime_error( "T temperature vector not found in atmosphere file" );
	}
	if (atm_profile->contains_vector("P")) {
		atm_profile->convert_property_units( "P", Units::fromString( "Pa" ) );
	} else {
		throw std::runtime_error( "P pressure vector not found in atmosphere file" );
	}
	if (atm_profile->contains_vector("RHO")) {
		atm_profile->convert_property_units( "RHO", Units::fromString( "kg/m3" ) );
	} else {
		throw std::runtime_error( "RHO density vector not found in atmosphere file" );
	}
	if (atm_profile->contains_scalar("Z0")) {
		atm_profile->convert_property_units( "Z0", Units::fromString( "m" ) );
	}

	// derived quantities
	if (atm_profile->contains_vector("WS")) {
		atm_profile->convert_property_units( "WS", Units::fromString("m/s") );
		atm_profile->copy_vector_property( "WS", "_WS_" );
	} else {
		atm_profile->calculate_wind_speed( "_WS_", "U", "V" );
	}
	if (atm_profile->contains_vector("WD")) {
		atm_profile->convert_property_units( "WD",
			Units::fromString("degrees clockwise from North") );
		atm_profile->copy_vector_property( "WD", "_WD_" );
	} else {
		atm_profile->calculate_wind_direction( "_WD_", "U", "V" );
	}

}


void NCPA::setup_ncpaprop_atmosphere( NCPA::Atmosphere2D *atm_profile_2d ) {

	// range units
	atm_profile_2d->convert_range_units( NCPA::Units::fromString( "m" ) );

	// altitude units
	atm_profile_2d->convert_altitude_units( Units::fromString( "m" ) );

	// static sound speed
	for (std::vector< NCPA::Atmosphere1D * >::iterator it = atm_profile_2d->first_profile();
		 it != atm_profile_2d->last_profile(); ++it) {
		if ( (*it)->contains_vector( "C0" ) ) {
			(*it)->convert_property_units( "C0", Units::fromString( "m/s" ) );
			(*it)->copy_vector_property( "C0", "_C0_" );
			// do we need to calculate temperature?
			if (!(*it)->contains_vector("T" )) {
				(*it)->calculate_temperature_from_sound_speed( "T", "_C0_",
					Units::fromString( "K" ) );
			}
		} else {
			if ( (*it)->contains_vector("P") && (*it)->contains_vector("RHO") ) {
				(*it)->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO",
					Units::fromString( "m/s" ) );
			} else if ( (*it)->contains_vector("T") ) {
				(*it)->calculate_sound_speed_from_temperature( "_C0_", "T",
					Units::fromString( "m/s" ) );
			} else if ( !((*it)->contains_vector("CEFF"))) {
				throw std::runtime_error( "Cannot calculate static sound speed: None of CEFF, C0, T, or (P and RHO) are specified in atmosphere.");
			}
		}

		// density
		if ( !(*it)->contains_vector( "RHO") ) {
			if ( (*it)->contains_vector("T") && (*it)->contains_vector("P") ) {
				(*it)->calculate_density_from_temperature_and_pressure(
					"RHO", "T", "P", Units::fromString( "kg/m3" ) );
			} else {
				throw std::runtime_error( "No RHO provided, and at least one of T and P is missing." );
			}
		}
	}

	// base parameters
	if (atm_profile_2d->contains_vector(0,"U")) {
		atm_profile_2d->convert_property_units( "U", Units::fromString( "m/s" ) );
	}
	if (atm_profile_2d->contains_vector(0,"V")) {
		atm_profile_2d->convert_property_units( "V", Units::fromString( "m/s" ) );
	}
	if (atm_profile_2d->contains_vector(0,"T")) {
		atm_profile_2d->convert_property_units( "T", Units::fromString( "K" ) );
	}
	if (atm_profile_2d->contains_vector(0,"P")) {
		atm_profile_2d->convert_property_units( "P", Units::fromString( "Pa" ) );
	}
	if (atm_profile_2d->contains_scalar(0,"Z0")) {
		atm_profile_2d->convert_property_units( "Z0", Units::fromString("m" ) );
	}
	if (atm_profile_2d->contains_vector(0,"RHO")) {
		atm_profile_2d->convert_property_units( "RHO", Units::fromString( "kg/m3" ) );
	}

	// derived quantities
	// wind speed
	if (atm_profile_2d->contains_vector(0,"WS")) {
		atm_profile_2d->convert_property_units("WS", Units::fromString("m/s") );
		atm_profile_2d->copy_vector_property( "WS", "_WS_" );
	} else if (atm_profile_2d->contains_vector(0,"U")
		&& atm_profile_2d->contains_vector(0,"V")) {
		atm_profile_2d->calculate_wind_speed( "_WS_", "U", "V" );
	}

	// wind direction
	if (atm_profile_2d->contains_vector(0,"WD")) {
		atm_profile_2d->convert_property_units("WD",
			Units::fromString("degrees clockwise from North") );
		atm_profile_2d->copy_vector_property( "WD", "_WD_" );
	} else if (atm_profile_2d->contains_vector(0,"U")
		&& atm_profile_2d->contains_vector(0,"V")) {
		atm_profile_2d->calculate_wind_direction( "_WD_", "U", "V" );
	}

}