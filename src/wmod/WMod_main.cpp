#include <cstdlib>
#include <cstring>

//#include "SolveWMod.h"
#include "modes.h"
#include "util.h"
#include "Atmosphere1D.h"
#include "wmod_parameters.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

/*
 * This is Jelle Assink's Normal Modes (wide angle)- moved to C++
 * Contains some of Roger Waxler's original design also.
 * @version 2.0
 * @date 2020-06
 * @authors Jelle Assink; Roger Waxler; Claus Hetzer; Doru Velea;
 * 
 * Changelog:
 * 20130326: DV added turnoff_WKB flag
 * 201305  : DV added use_attn_file option (to allow atten. coeff loaded from a text file)
 * 201306  : DV modified the get() functions that return strings
 * 202006  : CH retrofitted to use new atmospheric and input libraries
 */


//
// main
//
int main( int argc, char **argv ) {

  // object to process the options
  ParameterSet *param = new ParameterSet();
  try {
    configure_wmod_parameter_set( param );
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
  } catch (std::invalid_argument& e) {
    std::cout << "Parameter validation failed: " << std::endl
          << e.what() << std::endl;
  }

  // set up to measure the duration of this run
  time_t tm1 = time(NULL);

  // obtain the parameter values from the user's options
  // open the file
  string atmosfile = param->getString( "atmosfile" );
  string atmosheaderfile =  param->getString( "atmosheaderfile" );
  Atmosphere1D *atm_profile = new Atmosphere1D( atmosfile, atmosheaderfile );

  // get solver object 
  WModeSolver *a = new WModeSolver(param, atm_profile);
   
	//   					 
  // compute modes - main action happens here					 
  //
  try {
    a->solve();
  } catch (std::runtime_error &e) {
    std::cout << "Runtime error: " << e.what() << std::endl;
    return 1;
  }
  a->printParams();
  
  // save atm. profile if requested
  // if (param->getBool( "write_atm_profile" ) ) {
  //   ofstream ofs( "atm_profile.nm" );
  //   atm_profile->print_atmosphere( "Z", ofs );
  //   ofs.close();
  // } else {
  //   printf(" write_atm_profile flag : %d\n", param->getBool( "write_atm_profile" ));
  // }
  	  
  delete a;		  	 
  delete atm_profile;
  delete param;

	cout << endl << " ... main() is done." << endl;
	time_t tm2 = time(NULL);
  cout << "Run duration: " << difftime(tm2,tm1) << " seconds." << endl;	

	return 0;
} // end of main();


