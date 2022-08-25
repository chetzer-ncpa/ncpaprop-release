#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <string>
#include <vector>
#include <cfloat>
#include <fstream>
#include <stdexcept>

// timing
#include <chrono>

#include "petscksp.h"

#include "Atmosphere1D.h"
#include "parameterset.h"
#include "units.h"
#include "util.h"
//#include "epade_pe.h"
#include "EPadeSolver.h"
#include "epade_pe_parameters.h"

using namespace std;
using namespace NCPA;

int main( int argc, char **argv ) {
	
	using namespace std::chrono;
	
	std::string help = "ePade Parabolic Equation Solver";
	//PetscErrorCode ierr = PetscInitialize( &argc, &argv, (char *)0, help.c_str());CHKERRQ(ierr);
	PetscErrorCode ierr = PetscInitializeNoArguments();CHKERRQ(ierr);


	// object to process the options
	ParameterSet *param = new ParameterSet();
	configure_epade_pe_parameter_set( param );
	param->parseCommandLine( argc, argv );

	// check for help text
	if (param->wasFound( "help" ) || param->wasFound("h") ) {
		param->printUsage( cout );
		return 1;
	}

	// See if an options file was specified
	string paramFile = param->getString( "paramfile" );
	param->parseFile( paramFile );

	// parse command line again, to override file options
	param->parseCommandLine( argc, argv );

	// see if we want a parameter summary
	if (param->wasFound( "printparams" ) ) {
		param->printParameters();
	}

	// run parameter checks
	if (! param->validate() ) {
		cout << "Parameter validation failed:" << endl;
		param->printFailedTests( cout );
		return 0;
	}
	
	EPadeSolver *solver;
	try {
		// solver = new EPadeSolver( param );
		solver = new EPadeSolver();
		configure_epade_solver( solver, param );
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		solver->solve();
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		cout << "Elapsed time: " << time_span.count() << " seconds." << endl;
		delete solver;
	} catch (std::runtime_error& e) {
		std::cout << "ePape run failed with the following error:"
				  << endl << e.what() << endl;
	}

	/*
	cout << "Writing 1-D output to tloss_1d.pe" << endl;
	solver->output1DTL( "tloss_1d.pe" );
	cout << "Writing 2-D output to tloss_2d.pe" << endl;
	solver->output2DTL( "tloss_2d.pe" );
	*/

	// delete solver;
	delete param;

	ierr = PetscFinalize();

	return 1;
}