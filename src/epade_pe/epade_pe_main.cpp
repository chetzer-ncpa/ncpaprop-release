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

enum parameter_processing_results_t {
	CONTINUE_OK,
	EXIT_OK,
	EXIT_ERROR
};

void configure_solver( EPadeSolver *solver, ParameterSet *param );
parameter_processing_results_t process_parameters( ParameterSet *param, int argc, char **argv );


int main( int argc, char **argv ) {

	using namespace std::chrono;

	std::string help = "ePade Parabolic Equation Solver";
	//PetscErrorCode ierr = PetscInitialize( &argc, &argv, (char *)0, help.c_str());CHKERRQ(ierr);
	PetscErrorCode ierr = PetscInitializeNoArguments();CHKERRQ(ierr);


	// object to process the options
	ParameterSet *param = new ParameterSet();
	switch (process_parameters( param, argc, argv )) {
		case EXIT_OK:
			return 1;
		case EXIT_ERROR:
			return 0;
		default:
			break;
	}

	// validate parameters
	if (! param->validate() ) {
		cerr << "Parameter validation failed:" << endl;
		param->printFailedTests( cout );
		return 0;
	}
	
	EPadeSolver *solver;
	try {
		solver = new EPadeSolver( param );
		configure_solver( solver, param );
		solver->finalize();
		try {
			if (solver->ready()) {
				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				solver->solve();
				high_resolution_clock::time_point t2 = high_resolution_clock::now();
				duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
				cout << "Elapsed time: " << time_span.count() << " seconds." << endl;
			}
		} catch (std::runtime_error &e) {
			std::cerr << "Error in solver: " << e.what() << std::endl;
		}
		delete solver;
	} catch (std::runtime_error& e) {
		std::cout << "ePape run failed with the following error:"
				  << endl << e.what() << endl;
	}

	delete param;

	ierr = PetscFinalize();

	return 1;
}

parameter_processing_results_t process_parameters( ParameterSet *param, int argc, char **argv ) {
	configure_epade_pe_parameter_set( param );
	try {
		param->parseCommandLine( argc, argv );
	} catch (std::invalid_argument &e) {
		std::cout << "Argument parsing from command line failed:"
			<< std::endl << e.what() << std::endl;
		return EXIT_ERROR;
	}

	// check for help text
	if (param->wasFound( "help" ) || param->wasFound("h") ) {
		param->printUsage( cout );
		return EXIT_OK;
	}

	// See if an options file was specified
	string paramFile = param->getString( "paramfile" );
	try {
		param->parseFile( paramFile );
	} catch (std::invalid_argument &e) {
		std::cout << "Argument parsing from " << paramFile
			<< " failed: " << std::endl << e.what() << std::endl;
		return EXIT_ERROR;
	}

	// parse command line again, to override file options
	try {
		param->parseCommandLine( argc, argv );
	} catch (std::invalid_argument &e) {
		std::cout << "Argument parsing from command line failed:"
			<< std::endl << e.what() << std::endl;
		return EXIT_ERROR;
	}

	// see if we want a parameter summary
	if (param->wasFound( "printparams" ) ) {
		param->printParameters();
	}

	return CONTINUE_OK;
}


void configure_solver( EPadeSolver *solver, ParameterSet *param ) {
	solver->set_frequency( param->getFloat( "freq" ) );
	if (param->wasFound("maxheight_km")) {
		solver->set_max_height( param->getFloat( "maxheight_km" ), "km" );
	}
	if (param->wasFound("sourceheight_km")) {
		solver->set_source_height( param->getFloat( "sourceheight_km" ), "km" );
	}
	if (param->wasFound("receiverheight_km")) {
		solver->set_receiver_height( param->getFloat( "receiverheight_km" ), "km" );
	}
	if (param->wasFound("Nrng_steps")) {
		solver->set_requested_range_steps( param->getInteger( "Nrng_steps" ) );
	}
	if (param->wasFound("npade")) {
		solver->set_pade_order( param->getInteger( "npade" ) );
	}
	if (param->wasFound("dz_m")) {
		solver->set_requested_height_step( param->getFloat( "dz_m" ), "m" );
	}

	// starter
	std::string starter = param->getString("starter");
	if (starter == "self") {
		solver->set_starter( StarterType::SELF );
	} else if (starter == "gaussian") {
		solver->set_starter( StarterType::GAUSSIAN );
	} else if (starter == "user") {
		solver->set_starter( StarterType::USER, param->getString( "starterfile" ) );
	} else {
		throw std::runtime_error("Unrecognized starter type!");
	}

	// attenuation
	if (param->wasFound("lossless")) {
		solver->set_attenuation(AttenuationType::NONE);
	} else if (param->wasFound("attnfile") && param->getString("attnfile") != "") {
		solver->set_attenuation(AttenuationType::USER, param->getString("attnfile"));
	} else {
		solver->set_attenuation(AttenuationType::SUTHERLAND_BASS);
	}
}
