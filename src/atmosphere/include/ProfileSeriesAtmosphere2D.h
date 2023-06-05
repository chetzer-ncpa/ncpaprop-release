#ifndef NCPAPROP_PROFILESERIESATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_PROFILESERIESATMOSPHERE2D_H_INCLUDED

#include "Atmosphere2D.h"


namespace NCPA {
	class ProfileSeriesAtmosphere2D : public Atmosphere2D {

	public:
		ProfileSeriesAtmosphere2D();
		ProfileSeriesAtmosphere2D( const std::string &filename, const std::string &headerfilename = "" );
		~ProfileSeriesAtmosphere2D();

	};
}




#endif