#include "parameterset.h"
#include "modbb_parameters.h"
#include <string>
#include <iostream>

void NCPA::configure_modbb_parameter_set( NCPA::ParameterSet *ps ) {

	NCPA::ParameterTest *test = NULL;

	// general configuration
	// set up expected commands
	ps->setStrict( true );
	ps->setComments( "#%" );
	ps->setDelimiters( "=: " );

	// Add header instructions
	ps->addHeaderTextVerbatim( "----------------------------------------------------------------------------" );
	ps->addHeaderTextVerbatim( "|                             NCPA Infrasound                              |" );  	
	ps->addHeaderTextVerbatim( "|                         Normal Modes Broadband                           |" );
	ps->addHeaderTextVerbatim( "|    Based on either: Effective Sound Speed Approximation - see ModESS     |" );
	ps->addHeaderTextVerbatim( "|                     Wide_Angle High-Mach code - see WMod                 |" );  
	ps->addHeaderTextVerbatim( "----------------------------------------------------------------------------" );	
 	ps->addBlankHeaderLine();
	
 	ps->addHeaderText( "One of two algorithms can be used to perform pulse propagation.  The first is based on the Effective Sound Speed Approximation (as in Modess), the second is based on the the Wide_Angle High-Mach solution of the wave equation (as in WMod).  Modess is faster but is accurate for launch angles less than 30 degrees and low wind speeds.  WMod extends the validity to higher angles and high Mach numbers, but is slower.");
 	ps->addHeaderText( "To propagate a pulse, two steps must be completed:" );
 	ps->addHeaderText( "1. A dispersion file must be calculated using the option --dispersion ." );
 	ps->addHeaderText( "2. Pulse propagation is calculated for a selected source type: " );
 	ps->setHeaderIndent( 2 );
 	ps->addHeaderText( "--source impulse                           : Delta function" );
 	ps->addHeaderText( "--source pulse1                            : Built-in pulse type 1" );
 	ps->addHeaderText( "--source pulse2                            : Built-in pulse type 2" );
 	ps->addHeaderText( "--source spectrum --source_file <filename> : User-supplied spectrum file" );
 	ps->addHeaderText( "    Format: Freq   Re[ spec(f) ]  Im[ spec(f) ]" );
 	ps->addHeaderText( "--source waveform --source_file <filename> : User-supplied waveform file" );
 	ps->addHeaderText( "    Format: Time   Amplitude" );
 	ps->resetHeaderIndent();
	ps->addBlankHeaderLine();
	ps->addHeaderText("The options below can be specified in a colon-separated file \"modbb.param\" or at the command line. Command-line options override file options.");

	// Parameter descriptions
	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameter( new FlagParameter( "h" ) );
	ps->addParameterDescription( "Options Control", "--help", "Prints help test" );

	ps->addParameter( new StringParameter( "paramfile", "modbb.param") );
	ps->addParameterDescription( "Options Control", "--paramfile", "Parameter file name [modbb.param]" );

	ps->addParameter( new FlagParameter( "printparams" ) );
	ps->addParameterDescription( "Options Control", "--printparams", "Print parameter summary to screen" );

	// Modes of operation
	std::string modes_of_operation[ 2 ] = { "dispersion", "propagation" };
	for (unsigned int i = 0; i < 2; i++) {
		std::string tmpStr( modes_of_operation[ i ] );
		ps->addParameter( new NCPA::FlagParameter( tmpStr ) );
	}
	ps->addTest( new NCPA::RadioButtonTest( "operation_mode", 2, modes_of_operation ) );
	ps->addParameterDescription( "Modes of Operation", "--dispersion", "Calculate dispersion at multiple frequencies" );
	
	// required for dispersion calculation
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );

	ps->addParameter( new NCPA::StringParameter( "method", "modess" ) );
	test = ps->addTest( new NCPA::StringSetTest( "method" ) );
	test->addStringParameter( "modess" );
	test->addStringParameter( "wmod" );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "method", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--method", "Calculation method.  Currently supported: { modess, wmod } [required]" );

	ps->addParameter( new NCPA::StringParameter( "dispersion_file", "" ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "dispersion_file", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--dispersion_file", "File to write dispersion results to [required]" );

	ps->addParameter( new FlagParameter( "append_dispersion_file", true ) );
	  
	ps->addParameter( new NCPA::FloatParameter( "azimuth" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--azimuth", "Azimuth of propagation, in degrees CW from North [0,360) [required]" );
	
	ps->addParameter( new NCPA::FloatParameter( "f_min" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "f_min", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "f_min", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--f_min", "Minimum frequency in Hz [required]" );

	ps->addParameter( new NCPA::FloatParameter( "f_max" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "f_max", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "f_max", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--f_max", "Maximum frequency in Hz [required]" );

	ps->addParameter( new NCPA::FloatParameter( "f_step" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "f_step", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "f_step", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--f_step", "Frequency interval in Hz [required]" );

	ps->addParameter( new NCPA::StringParameter( "atmosfile" ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "atmosfile", "dispersion" ) );
	ps->addParameterDescription( "Modes of Operation", "--atmosfile", "Atmospheric profile filename [required]" );

	ps->addParameter( new NCPA::StringParameter( "atmosheaderfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--atmosheaderfile", "External header file, overrides internal header [None]" );

	ps->addParameter( new NCPA::FloatParameter( "maxheight_km", 150.0 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxheight_km", 0.1 ));
	ps->addParameterDescription( "Modes of Operation", "--maxheight_km", "Maximum height of analysis in km [150.0]" );


	ps->addParameter( new NCPA::FloatParameter( "zground_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--zground_km", "Ground height [take from Z0 parameter in atmosfile, or 0.0]" );


	ps->addParameter( new NCPA::IntegerParameter( "Nz_grid", 20000 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "Nz_grid", 10 ) );
	ps->addParameterDescription( "Modes of Operation", "--Nz_grid", "Number of vertical grid points to use [20000]" );

	ps->addParameter( new NCPA::FloatParameter( "sourceheight_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--sourceheight_km", "Source height in km [0.0]" );


	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--receiverheight_km", "Receiver height in km [0.0]" );

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km", 1000.0 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxrange_km", 0.01 ) );
	ps->addParameterDescription( "Modes of Operation", "--maxrange_km", "Maximum range in km to use for modeling [1000.0]" );

	ps->addParameter( new NCPA::IntegerParameter( "Nrng_steps", 1000 ) );
	ps->addTest( new NCPA::IntegerGreaterThanTest( "Nrng_steps", 0 ) );
	ps->addParameterDescription( "Modes of Operation", "--Nrng_steps", "Number of range steps to use [1000]" );

	ps->addParameter( new NCPA::StringParameter( "ground_impedence_model", "rigid" ) );
	test = ps->addTest( new NCPA::StringSetTest( "ground_impedence_model" ) );
	test->addStringParameter( "rigid" );
	ps->addParameterDescription( "Modes of Operation", "--ground_impedence_model", "Impedence model to use.  Currently only \"rigid\" is supported. [rigid]" );

	ps->addParameter( new NCPA::FlagParameter( "Lamb_wave_BC" ) );
	ps->addParameterDescription( "Modes of Operation", "--Lamb_wave_BC", "Use admittance = -1/2*dln(rho)/dz" );

	ps->addParameter( new NCPA::FlagParameter( "write_atm_profile" ) );
	ps->addParameterDescription( "Modes of Operation", "--write_atm_profile", "Output atmospheric profile to atm_profile.nm" );

	ps->addParameter( new NCPA::StringParameter( "use_attn_file", "" ) );
	ps->addParameterDescription( "Modes of Operation", "--use_attn_file", "File name containing attenuation, to override default Sutherland/Bass [n/a]. Columns are #n# Height(km) Attenuation(np/m)" );

	ps->addParameter( new NCPA::FlagParameter( "use_zero_attenuation" ) );
	ps->addParameterDescription( "Modes of Operation", "--use_zero_attenuation", "Set attenuation to zero." );

	ps->addParameter( new NCPA::FlagParameter( "wvnum_filter" ) );
	ps->addParameterDescription( "Modes of Operation", "--wvnum_filter", "Use wavenumber filter by phase speed.  Requires --c_min and --c_max" );

	ps->addParameter( new NCPA::FloatParameter( "c_min" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "c_min", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "c_min", "wvnum_filter" ) );
	ps->addParameter( new NCPA::FloatParameter( "c_max" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "c_max", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "c_max", "wvnum_filter" ) );
	ps->setParameterIndent( 3 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--c_min", "Minimum phase speed to keep" );
	ps->addParameterDescription( "Modes of Operation", "--c_max", "Maximum phase speed to keep" );
	ps->resetParameterIndent();

	// options for propagation of calculated dispersion
	ps->addParameterDescription( "Modes of Operation", "--propagation", "Calculate propagated waveform using dispersion results" );
	
	// required for dispersion calculation
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );

	ps->addParameter( new NCPA::StringParameter( "input_dispersion_file", "" ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "input_dispersion_file", "propagation" ) );
	ps->addParameterDescription( "Modes of Operation", "--input_dispersion_file", "File to read dispersion results from [required]" );

	ps->addParameter( new NCPA::StringParameter( "output_waveform_file", "" ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "output_waveform_file", "propagation" ) );
	ps->addParameterDescription( "Modes of Operation", "--output_waveform_file", "File to write resultant waveform to [required]" );

	ps->addParameter( new NCPA::StringParameter( "source", "impulse" ) );
	test = ps->addTest( new NCPA::StringSetTest( "source" ) );
	test->addStringParameter( "impulse" );
	test->addStringParameter( "pulse1" );
	test->addStringParameter( "pulse2" );
	test->addStringParameter( "waveform" );
	test->addStringParameter( "spectrum" );
	ps->addParameterDescription( "Modes of Operation", "--source", "Source type.  Options include: {impulse,pulse1,pulse2,spectrum,waveform} [impulse]" );

	ps->setParameterIndent( 3 * DEFAULT_PARAMETER_INDENT );
	ps->addParameter( new NCPA::StringParameter( "source_file", "" ) );
	ps->addParameterDescription( "Modes of Operation", "--source_file", "File containing the source spectrum or waveform, if applicable" );
	ps->addParameter( new NCPA::FloatParameter( "f_center", -1.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--f_center", "Center frequency for pulse1 or pulse2 options.  Must be <= f_max/5 [f_max/5]");
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );

	ps->addParameter( new NCPA::IntegerParameter( "nfft", 0 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "nfft", 0 ) );
	ps->addParameterDescription( "Modes of Operation", "--nfft", "Number of FFT points [4*f_max/f_step]");
	
	ps->addParameter( new NCPA::FloatParameter( "max_celerity", 340.0 ) );
	ps->addTest( new NCPA::IntegerGreaterThanTest( "max_celerity", 0 ) );
	ps->addParameterDescription( "Modes of Operation", "--max_celerity", "Maximum celerity for calculation [340.0]");
	
	ps->addParameter( new NCPA::StringParameter( "receiver", "single" ) );
	test = ps->addTest( new NCPA::StringSetTest( "receiver" ) );
	test->addStringParameter( "single" );
	test->addStringParameter( "multiple" );
	ps->addParameterDescription( "Modes of Operation", "--receiver", "Receiver type {single,multiple} [single]" );

	ps->setParameterIndent( 3 * DEFAULT_PARAMETER_INDENT );
	ps->addParameter( new NCPA::FloatParameter( "range_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--range_km", "Propagation range in km for a single receiver");
	ps->addParameter( new NCPA::FloatParameter( "start_range_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--start_range_km", "Starting propagation range in km for multiple receivers");
	ps->addParameter( new NCPA::FloatParameter( "end_range_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--end_range_km", "Ending propagation range in km for multiple receivers");
	ps->addParameter( new NCPA::FloatParameter( "range_step_km", 0.0 ) );
	ps->addParameterDescription( "Modes of Operation", "--range_step_km", "Propagation range step in km for multiple receivers");
	ps->resetParameterIndent();

	// automatically broadband mode
	ps->addParameter( new NCPA::FlagParameter( "broadband", true ) );

	
	// dummy parameters, required but not used.  Need to handle this better in ESSModeSolver()
	// and others, but apparently can't do it with exceptions b/c clang.  Second-order problem
	// to solve.
	ps->addParameter( new NCPA::StringParameter( "modal_starter_file", "" ) );
	ps->addParameter( new NCPA::FlagParameter( "write_2d_tloss" ) );
	ps->addParameter( new NCPA::FlagParameter( "write_phase_speeds" ) );
	ps->addParameter( new NCPA::FlagParameter( "write_speeds" ) );
	ps->addParameter( new NCPA::FlagParameter( "write_modes" ) );
	ps->addParameter( new NCPA::FlagParameter( "write_lossless" ) );
	ps->addParameter( new NCPA::FlagParameter( "multiprop" ) );
	ps->addParameter( new NCPA::FlagParameter( "turnoff_WKB" ) );
	ps->addParameter( new NCPA::StringParameter( "filetag", "" ) );
	

	// Footer with file formats and sample commands
	ps->addBlankFooterLine();
	ps->addFooterText("OUTPUT Files:  Format description (column order):");
	ps->addFooterTextVerbatim("  <dispersion file>        Contains one line per frequency with entries:");
	ps->addFooterTextVerbatim("                           freq, (# of modes), rho(z_src),");
	ps->addFooterTextVerbatim("                           followed for each mode 'i' by quadruples:");
	ps->addFooterTextVerbatim("                           real(k(i)), imag(k(i)), Mode(i)(z_src), Mode(i)(z_rcv)");
	ps->addFooterTextVerbatim("  <waveform file>          r[km]  t[s]  P" );
	ps->addFooterText("Examples (run from 'samples' directory):");
	ps->setFooterIndent( 4 );
	ps->setFooterHangingIndent( 4 );
	ps->setCommandMode( true );
	ps->addFooterText("../bin/ModBB --dispersion --dispersion_file myDispersionFile.dat --atmosfile NCPA_canonical_profile_zuvwtdp.dat --azimuth 90 --f_min 0.001953125 --f_step 0.001953125 --f_max 0.5 --method modess" );
	ps->addFooterTextVerbatim( " --or--" );
	ps->addFooterText("../bin/ModBB --dispersion --dispersion_file myDispersionFile.dat --atmosfile NCPA_canonical_profile_zuvwtdp.dat --azimuth 90 --f_min 0.001953125 --f_step 0.001953125 --f_max 0.5 --method wmod" );
	ps->addBlankFooterLine();
	ps->addFooterTextVerbatim( " --then--" );
	ps->addFooterText("../bin/ModBB --propagation --input_dispersion_file myDispersionFile.dat --range_km 240 --output_waveform_file mywavf.dat --receiver single --source impulse" );
	ps->addBlankFooterLine();
	ps->addFooterText("../bin/ModBB --propagation --input_dispersion_file myDispersionFile.dat --output_waveform_file mywavf.dat --receiver multiple --start_range_km 240 --end_range_km 300 --range_step_km 20 --source waveform --source_file source_waveform_input_example.dat" );
	ps->setFooterHangingIndent( 0 );
	ps->setCommandMode( false );
	ps->resetFooterIndent();
}
