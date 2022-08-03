#ifndef NCPAPROP_ATMOSPHERE_INTERFACE_H_INCLUDED
#define NCPAPROP_ATMOSPHERE_INTERFACE_H_INCLUDED

#include "Atmosphere1D.h"
#include "Atmosphere2D.h"
#include "Atmosphere3D.h"

// functions to interface between the atmosphere library and ncpaprop
// modules
namespace NCPA {

	void setup_ncpaprop_atmosphere( NCPA::Atmosphere1D *profile );
	void setup_ncpaprop_atmosphere( NCPA::Atmosphere2D *profile );
	void setup_ncpaprop_atmosphere( NCPA::Atmosphere3D *profile );

}

#endif