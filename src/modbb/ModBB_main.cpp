#include <complex>
#include <string>
#include <ctime>
#include "modbb_parameters.h"
#include "ModeSolver.h"
#include "ESSModeSolver.h"
#include "WModeSolver.h"
#include "parameterset.h"
#include "Atmosphere1D.h"
#include "BroadbandPropagator.h"
#include "ModalBroadbandPropagator.h"


/*
 * This is the Normal Modes Broadband code 
 * based on either the Effective Sound Speed Approximation (Modess) or
 * on the Wide-Angle High-Mach modal code (WMod)
 * @authors Doru Velea; Jelle Assink; Roger Waxler; Claus Hetzer; 
 */

using namespace NCPA;
using namespace std;

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000


//
// main
//
int main( int argc, char **argv ) {

  // object to process the options
  ParameterSet *param = new ParameterSet();
  configure_modbb_parameter_set( param );
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

  // set the timer to measure running time
  std::time_t tm1 = std::time(NULL);

  //
  // logic to handle either computing the dispersion and modal values 
  // or to propagate a pulse
  //
  if (param->getBool( "dispersion" )) { //compute and write dispersion
      string atmosfile = param->getString( "atmosfile" );
      string atmosheaderfile =  param->getString( "atmosheaderfile" );
      Atmosphere1D *atm_profile = new Atmosphere1D( atmosfile, atmosheaderfile );

      // initialize dispersion file
      std::string dispersion_file = param->getString( "dispersion_file" );
      if (!dispersion_file.empty()) {
        std::ofstream dfs( dispersion_file, std::ios_base::trunc );
        dfs.close();
      }

      ModeSolver *solver;
      
      if (param->getString("method") == "modess") {
        solver = new ESSModeSolver( param, atm_profile );
      } else {
        solver = new WModeSolver( param, atm_profile );
      }
      solver->solve();
      delete solver;
      
      // save atm. profile if requested
      // if (param->getBool( "write_atm_profile" ) ) {
      //   ofstream ofs( "atm_profile.nm" );
      //   atm_profile->print_atmosphere( "Z", ofs );
      //   ofs.close();
      // }

      delete atm_profile;
      cout << "Finished computing and saving dispersion data" << endl;	
								
  } // end computing dispersion data

  else if (param->getBool( "propagation" )) { //pulse propagation source-to-receiver

    BroadbandPropagator *prop = new ModalBroadbandPropagator( param );
    prop->calculate_waveform();
    delete prop;

  }

  delete param; // delete object processing options

  std::time_t tm2 = std::time(NULL);
  cout << "\nRun duration: " << difftime(tm2,tm1) << " seconds." << endl;
  cout << " ... main() broadband version is done." << endl;

  return 0;
} // end of Broadband main();