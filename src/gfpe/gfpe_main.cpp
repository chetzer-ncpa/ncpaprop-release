#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <complex>
#include "matrix.h"
#include "gfpe_parameters.h"
#include "GFPESolver.h"

// timing
#include <chrono>

#ifndef PI
#define PI 3.14159265358979323846
#endif


using namespace std;
using namespace NCPA;


int main( int argc, char **argv ) {

	using namespace std::chrono;

	// Process input options
	ParameterSet *paramset = new ParameterSet();
	configure_gfpe_parameter_set( paramset );
	paramset->parseCommandLine( argc, argv );

	// check for help text
	if (paramset->wasFound( "help" ) || paramset->wasFound("h") ) {
		paramset->printUsage( cout );
		return 1;
	}

	// See if an options file was specified
	string paramFile = paramset->getString( "paramfile" );
	paramset->parseFile( paramFile );

	// parse command line again, to override file options
	paramset->parseCommandLine( argc, argv );

	// see if we want a parameter summary
	if (paramset->wasFound( "printparams" ) ) {
		paramset->printParameters();
	}

	// run parameter checks
	if (! paramset->validate() ) {
		cout << "Parameter validation failed:" << endl;
		paramset->printFailedTests( cout );
		return 0;
	}


	GFPESolver *solver;
	try {
		// run the solver
		solver = new GFPESolver( paramset );
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		solver->solve();
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		cout << "Elapsed time: " << time_span.count() << " seconds." << endl;

		// output 2-D transmission loss to file
		if (paramset->wasFound( "write_2d_tloss" ) ) {
			solver->truncate2DFile();
			cout << "Writing 2-D signal field..." << endl;
			for (size_t i = 0; i < solver->runs(); i++) {
				solver->output2DTL( i, true );
			}
		}

		/*
		Example of how to retrieve 1D signal level programmatically:

		cout << endl << "Retrieving results at z~=1 via get1DTL():" << endl;
		vector<double> r;
		vector<complex<double>> tl;
		double z_desired = 1.0;
		for (size_t i = 0; i < solver->runs(); i++) {
			solver->get1DTL( i, z_desired, r, tl );
			for (size_t j = 0; j < r.size(); j++) {
				cout << "r == " << r[j] << ": " << tl[j] << endl;
			}
		}


		Example of how to retrieve 2D signal level programmatically:

		cout << endl << "Retrieving results via get2DTL():" << endl;
		vector<double> z;
		NCPA::DenseMatrix<complex<double>> *tl_mat;
		for (size_t i = 0; i < solver->runs(); i++) {
			solver->get2DTL( i, r, z, tl_mat );
			for (size_t j = 0; j < r.size(); j++) {
				for (size_t k = 0; k < z.size(); k++) {
					cout << "r = " << r[j] << ", z = " << z[k]
						 << ": " << tl_mat->at( 0, 0 ) << endl;
				}
			}
			delete tl_mat;
		}
		*/

		// clean up
		delete paramset;
		delete solver;

	} catch (runtime_error &e) {
		cout << "GFPE run failed with the following error:"
			 << endl << e.what() << endl;

	}




	return 1;



}