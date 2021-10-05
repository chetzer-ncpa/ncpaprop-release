#ifndef NCPAPROP_STRATIFIEDATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_STRATIFIEDATMOSPHERE2D_H_INCLUDED

#include "Atmosphere2D.h"


namespace NCPA {
	class StratifiedAtmosphere2D : public Atmosphere2D {

	public:
		StratifiedAtmosphere2D( const Atmosphere1D *atm );
		StratifiedAtmosphere2D( const std::string &filename,
			const std::string &headerfilename = "" );
		//StratifiedAtmosphere2D( const StratifiedAtmosphere2D &atm );
		~StratifiedAtmosphere2D();

		virtual double get_interpolated_ground_elevation( double range );
		virtual double get_interpolated_ground_elevation_first_derivative( double range );
		virtual double get_interpolated_ground_elevation_second_derivative( double range );
	};
}




#endif