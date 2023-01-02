#include "StratifiedAtmosphere2D.h"
#include "Atmosphere1D.h"
#include "units.h"


NCPA::StratifiedAtmosphere2D::StratifiedAtmosphere2D( const Atmosphere1D *atm ) : Atmosphere2D() {
	Atmosphere1D *tempatm = new Atmosphere1D( *atm );
	set_insert_range_units( NCPA::Units::fromString( "km" ) );
	insert_profile( tempatm, 0.0 );
	sorted_ = true;
	delete tempatm;
}

NCPA::StratifiedAtmosphere2D::StratifiedAtmosphere2D( const std::string &filename,
		const std::string &headerfilename ) : Atmosphere2D() {
	Atmosphere1D *tempatm = new Atmosphere1D( filename, headerfilename );
	set_insert_range_units( NCPA::Units::fromString( "km" ) );
	insert_profile( tempatm, 0.0 );
	sorted_ = true;
	delete tempatm;
}

NCPA::StratifiedAtmosphere2D::~StratifiedAtmosphere2D() { }

// double NCPA::StratifiedAtmosphere2D::get_interpolated_ground_elevation( double range ) {
// 	return profiles_[0]->get("Z0");
// }

// double NCPA::StratifiedAtmosphere2D::get_interpolated_ground_elevation_first_derivative( double range ) {
// 	return 0;
// }

// double NCPA::StratifiedAtmosphere2D::get_interpolated_ground_elevation_second_derivative( double range ) {
// 	return 0;
// }

//NCPA::StratifiedAtmosphere2D::StratifiedAtmosphere2D( const StratifiedAtmosphere2D &atm ) : Atmosphere2D( atm ) { }
