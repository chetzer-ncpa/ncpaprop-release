#ifndef NCPAPROP_EPADE_PE_PARAMETERS_H_DEFINED
#define NCPAPROP_EPADE_PE_PARAMETERS_H_DEFINED

#include "parameterset.h"
#include "EPadeSolver.h"
#include <iostream>

namespace NCPA {
	void configure_epade_pe_parameter_set( NCPA::ParameterSet *params );
	void configure_epade_solver( NCPA::EPadeSolver *solver,
		NCPA::ParameterSet *params );
}

#endif