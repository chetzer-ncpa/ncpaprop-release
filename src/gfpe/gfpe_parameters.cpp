#include "gfpe_parameters.h"
#include "parameterset.h"
#include "units.h"
#include <string>
#include <iostream>

void NCPA::configure_gfpe_parameter_set( NCPA::ParameterSet *ps ) {

	// general configuration
	ps->setStrict( true );
	ps->setComments( "#%" );
	ps->setDelimiters( "=: " );

	// Add header instructions
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addHeaderTextVerbatim("|                             NCPA Infrasound                              |");
	ps->addHeaderTextVerbatim("|                      Green's Function Parabolic Equation                 |");
	ps->addHeaderTextVerbatim("|          Multi-Frequency - Effective Sound Speed Approximation           |");
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addBlankHeaderLine();
	ps->addHeaderText("The program computes both the 1D (at the ground) and 2D transfer function and saves the data to 2 files:" );
	ps->setHeaderIndent( 4 );
	ps->addHeaderText("file tloss_1d.pe - transfer function at the ground (default)" );
	ps->addHeaderText("file tloss_2d.pe - full 2-D transfer function field (on request)" );
	ps->resetHeaderIndent();
	ps->addBlankHeaderLine();
	ps->addHeaderText("The NCPA GFPE can be used in one of two modes, one-run mode or suite mode.");
	ps->addBlankHeaderLine();
	ps->addHeaderText("One-run mode is a standard PE calculation where the 2-D signal level field in (r,z) coordinates is calculated once per frequency-azimuth pair, with or without random turbulence, and the full field retained for later retrieval (using the get1DTL() and get2DTL() methods) and/or file output (using the output1DTL() and output2DTL() methods).  Note that, as in the ePape Pade PE module, 1-D signal level file output happens automatically to the default filename tloss_1d.pe.  A unique tag can be prepended to the filename using the --filetag parameter.");
	ps->addBlankHeaderLine();
	ps->addHeaderText("Suite mode runs a specified number of identical (apart from random turbulence) calculations for each frequency-azimuth pair, and retains signal levels only at (r,z) pairs specified by the user.  If no r or z coordinates are specified, they default to ground level at the maximum range of the calculation.  These are retrieved and output in the same way as in one-run mode, but with output at specified coordinates only.  Suite mode is activated by specifying a value greater than 1 for the --n_suite parameter.  Suite mode can be used without turbulence, but will return identical signal levels for each run.");
	ps->addBlankHeaderLine();
	ps->addHeaderText("Examples of how to retrieve, display, and output signal leves are provided in the commented areas of gfpe.cpp.");
	ps->addHeaderText("The options below can be specified in a colon-separated file \"gfpe.param\" or at the command line. Command-line options override file options.");

	// Parameter descriptions
	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameter( new FlagParameter( "h" ) );
	ps->addParameterDescription( "Options Control", "--help", "Prints help test" );

	ps->addParameter( new StringParameter( "paramfile", "gfpe.param") );
	ps->addParameterDescription( "Options Control", "--paramfile", "Parameter file name [gfpe.param]" );

	ps->addParameter( new FlagParameter( "printparams" ) );
	ps->addParameterDescription( "Options Control", "--printparams", "Print parameter summary to screen" );

	// Atmosphere
	std::string atmosphere_types[ 2 ] = { "atmosfile", "atmosfile2d" };
	ps->addParameter( new NCPA::StringParameter( atmosphere_types[ 0 ] ) );
	ps->addParameter( new NCPA::StringParameter( atmosphere_types[ 1 ] ) );
	ps->addTest( new NCPA::RadioButtonTest( "atmosphere_type", 2, atmosphere_types ) );
	ps->addParameterDescription( "Atmosphere", "--atmosfile", "1-D atmospheric profile filename" );
	ps->addParameterDescription( "Atmosphere", "--atmosfile2d", "2-D atmospheric summary filename (see NCPAProp manual).  Note that non-stratified atmospheres are not yet supported." );

	ps->addParameter( new NCPA::FloatParameter( "humidity", 45 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "humidity", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "humidity", 100.0 ) );
	ps->addParameterDescription( "Atmosphere", "--humidity", "Relative humidity (percent) [45]" );

	// required parameters
	// Frequency
	ps->addParameter( new NCPA::FloatParameter( "freq" ) );
	ps->addTest( new NCPA::RequiredTest( "freq" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "freq", 0.0 ) );
	ps->addParameterDescription( "Required Parameters", "--freq", "Single frequency (Hz)" );

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km" ) );
	ps->addTest( new NCPA::RequiredTest( "maxrange_km" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxrange_km", 0.001 ) );
	ps->addParameterDescription( "Required Parameters", "--maxrange_km", "Maximum range in km to use for modeling" );



	// Modes of operation
	std::string modes_of_operation[ 2 ] = { "singleprop", "multiprop" };
	for (unsigned int i = 0; i < 2; i++) {
		std::string tmpStr( modes_of_operation[ i ] );
		ps->addParameter( new NCPA::FlagParameter( tmpStr ) );
	}
	ps->addTest( new NCPA::RadioButtonTest( "operation_mode", 2, modes_of_operation ) );
	ps->addParameterDescription( "Modes of Operation", "--singleprop", "Single azimuth propagation.  Requires --azimuth" );

	// for single propagation, must specify azimuth
	ps->addParameter( new NCPA::FloatParameter( "azimuth" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth", -360.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth", "singleprop" ) );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth", "Propagation azimuth ( degrees clockwise from north, [0,360) )" );
	ps->resetParameterIndent();

	// for multiple propagation, must specify azimuth start, end, and step
	ps->addParameterDescription( "Modes of Operation", "--multiprop", "Multiple azimuth propagation.  Requires --azimuth_start, --azimuth_end, and --azimuth_step, disables --atmosfile2d.  If azimuth_start > azimuth_end (e.g. to wrap around North), azimuth_start will be adjusted to its negative-azimuth equivalent." );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_start" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_start", -360.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_start", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_start", "multiprop" ) );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_end" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_end", -360.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_end", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_end", "multiprop" ) );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_step" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_step", -360.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_step", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_step", "multiprop" ) );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_start", "Starting azimuth, in degrees CW from North [0,360)" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_end", "Ending azimuth, in degrees CW from North [0,360)" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_step", "Azimuth step, in degrees CW from North [0,360)" );
	ps->resetParameterIndent();




	// optional parameters
	ps->addParameter( new NCPA::IntegerParameter( "n_suite", 1 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--n_suite", "Number of individual runs to compute [1]" );

	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameter( new NCPA::StringParameter( "r_suite_km" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--r_suite_km", "For n_suite > 1, comma-separated list of ranges to output, in km [end range only]" );
	ps->addParameter( new NCPA::StringParameter( "z_suite_km" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--z_suite_km", "For n_suite > 1, comma-separated list of heights to output, in km [ground height only]" );
	ps->resetParameterIndent();

	ps->addParameter( new NCPA::StringParameter( "atmosheaderfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--atmosheaderfile", "External header file, overrides internal header [None]" );

	ps->addParameter( new NCPA::FloatParameter( "sourceheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--sourceheight_km", "Source height in km [ground]" );

	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--receiverheight_km", "Receiver height in km [ground]" );

	ps->addParameter( new NCPA::FloatParameter( "humidity", -1.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--humidity", "Relative humidity.  Overrides any 'H' field in atmospheric file.  Humidity is required if density is not supplied." );

	ps->addParameter( new NCPA::IntegerParameter( "n_turbulence", 20 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--n_turbulence", "Number of random turbulence phases to compute [20]" );

	ps->addParameter( new NCPA::StringParameter( "turbulence_file" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_file", "File containing 2*n_turbulence numbers in [0,1) range [randomize]" );

//	ps->addParameter( new NCPA::FloatParameter( "turbulence_ref_temp", 293.0 ) );
//	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_ref_temp", "Reference temperature, in K, for turbulence calculation [293]" );

	ps->addParameter( new NCPA::FloatParameter( "turbulence_scale_m", 100.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_scale_m", "Turbulence scale (m) [100]" );

	ps->addParameter( new NCPA::FloatParameter( "turbulence_t_factor", 1.0e-10 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_t_factor", "Temperature factor for turbulence spectrum (m^-(2/3)) [1.0e-10]" );

	ps->addParameter( new NCPA::FloatParameter( "turbulence_v_factor", 1.0e-8 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_v_factor", "Wind velocity factor for turbulence spectrum (m^-(2/3)) [1.0e-8]" );

	ps->addParameter( new NCPA::FloatParameter( "k_min", 0.1 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--k_min", "Wavenumber filter minimum [0.1]" );

	ps->addParameter( new NCPA::FloatParameter( "k_max", 0.9 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--k_max", "Wavenumber filter maximum [0.9]" );

	ps->addParameter( new NCPA::FloatParameter( "ground_impedence", 200000.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence", "Flow resistivity for ground impedence calculation [200000]" );

	ps->addParameter( new NCPA::IntegerParameter( "nz", -1 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--nz", "Minimum number of vertical points [internal]" );

	ps->addParameter( new NCPA::FloatParameter( "surface_layer_m", -1.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--surface_layer_m", "Thickness of surface layer (m) [100 wavelengths]" );

	ps->addParameter( new NCPA::FloatParameter( "boundary_layer_m", -1.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--boundary_layer_m", "Thickness of boundary layer (m) [50 wavelengths]" );

	ps->addParameter( new NCPA::FloatParameter( "dr_m", -1.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--dr_m", "Horizontal step size (m) [10 wavelengths]" );

	ps->addParameter( new NCPA::StringParameter( "filetag", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--filetag", "Tag prepended to output filenames [n/a]" );


	// flags
	ps->addParameter( new NCPA::FlagParameter( "write_2d_tloss" ) );
	ps->addParameterDescription( "Flags", "--write_2d_tloss", "Output 2-D transfer function to tloss_2d.pe" );
	ps->addParameter( new NCPA::FlagParameter( "no_turbulence" ) );
	ps->addParameterDescription( "Flags", "--no_turbulence", "Disable turbulence perturbation" );
	ps->addParameter( new NCPA::FlagParameter( "no_delay" ) );
	ps->addParameterDescription( "Flags", "--no_delay", "Do not apply propagation phase delay" );
	ps->addParameter( new NCPA::FlagParameter( "quiet" ) );
	ps->addParameterDescription( "Flags", "--quiet", "Suppress output" );

	ps->addBlankFooterLine();
	ps->addFooterText("OUTPUT Files:  Format description (column order):");
	ps->addFooterTextVerbatim("  tloss_1d.pe:                 r (km), az (deg), TF (real), TF (imag)");
	// ps->addFooterTextVerbatim("  tloss_multiplot.pe:          r (km), az (deg), TF (real), TF (imag)");
	ps->addFooterTextVerbatim("  tloss_2d.pe:                 az, r, z, TF (real), TF (imag)");
	// ps->addFooterTextVerbatim("  atm_profile.pe:              z,u,v,w,t,d,p,c,c_eff");
	// ps->addFooterTextVerbatim("  starter.pe:                  z, starter (real), starter(imag)" );
	// ps->addFooterTextVerbatim("  topography.pe:               az, r, z0" );
	ps->addBlankFooterLine();
	ps->addFooterText("Examples (run from samples directory):");
	ps->setFooterIndent( 4 );
	ps->setFooterHangingIndent( 4 );
	ps->setCommandMode( true );
	ps->addBlankFooterLine();

	ps->addFooterText("../bin/gfpe --atmosfile synthetic.atm --singleprop --azimuth 90 --maxrange_km 0.5 --sourceheight_km 0.001 --receiverheight_km 0.001 --write_2d_tloss --freq 800" );
	ps->addBlankFooterLine();
	ps->addFooterText("../bin/gfpe --atmosfile synthetic.atm --multiprop --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 --freq 400 --maxrange_km 0.5 --sourceheight_km 0.001 --receiverheight_km 0.001" );
	ps->addBlankFooterLine();
	ps->addFooterText("../bin/gfpe --atmosfile synthetic.atm --singleprop --azimuth 90 --maxrange_km 1.0 --freq 500 --n_suite 50");
	ps->addBlankFooterLine();
	
	ps->setFooterHangingIndent( 0 );
	ps->setCommandMode( false );
	ps->resetFooterIndent();

}
