#include "parameterset.h"
#include "wmod_parameters.h"
#include <string>


void NCPA::configure_wmod_parameter_set( NCPA::ParameterSet *ps ) {

	NCPA::ParameterTest *test = NULL;

	// general configuration
	// set up expected commands
	ps->setStrict( true );
	ps->setComments( "#%" );
	ps->setDelimiters( "=: " );

	// Add header instructions
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addHeaderTextVerbatim("|                             NCPA Infrasound                              |");
	ps->addHeaderTextVerbatim("|                Wide Angle High-Mach-Number Normal Modes                  |");
	ps->addHeaderTextVerbatim("|                            Single Frequency                              |");
	ps->addHeaderTextVerbatim("|                    Attenuation added perturbatively                      |");
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addBlankHeaderLine();
	ps->addHeaderText("This program computes the 1D transmission loss (TL) using a wide-angle high-mach-number normal mode technique.  The computed loss is output in one or more files as controlled by user-supplied flags.  Default filenames and their contents are detailed below." );
	ps->addBlankHeaderLine();
	ps->addHeaderText("The options below can be specified in a colon-separated file \"wmod.param\" or at the command line. Command-line options override file options.");

	// Parameter descriptions
	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameter( new FlagParameter( "h" ) );
	ps->addParameterDescription( "Options Control", "-h, --help", "Prints help test" );
	ps->addParameter( new StringParameter( "paramfile", "modess.param") );
	ps->addParameterDescription( "Options Control", "--paramfile", "Parameter file name [modess.param]" );
	ps->addParameter( new FlagParameter( "printparams" ) );
	ps->addParameterDescription( "Options Control", "--printparams", "Print parameter summary to screen" );

	// Required parameters
	ps->addParameter( new NCPA::StringParameter( "atmosfile" ) );
	ps->addTest( new NCPA::RequiredTest( "atmosfile" ) );
	ps->addParameterDescription( "Required Parameters", "--atmosfile", "Atmospheric profile filename" );

	ps->addParameter( new NCPA::FloatParameter( "freq" ) );
	ps->addTest( new NCPA::RequiredTest( "freq" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "freq", 0.0 ) );
	ps->addParameterDescription( "Required Parameters", "--freq", "Frequency of analysis (Hz)" );

	// optional parameters
	ps->addParameter( new NCPA::StringParameter( "atmosheaderfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--atmosheaderfile", "External header file, overrides internal header [None]" );

	ps->addParameter( new NCPA::FloatParameter( "maxheight_km", 150.0 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxheight_km", 0.1 ));
	ps->addParameterDescription( "Optional Parameters [default]", "--maxheight_km", "Maximum height of analysis in km [150.0]" );


	ps->addParameter( new NCPA::FloatParameter( "zground_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--zground_km", "Ground height [take from Z0 parameter in atmosfile, or 0.0]" );


	ps->addParameter( new NCPA::IntegerParameter( "Nz_grid", 20000 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "Nz_grid", 10 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nz_grid", "Number of vertical grid points to use [20000]" );

	ps->addParameter( new NCPA::FloatParameter( "sourceheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--sourceheight_km", "Source height in km [0.0]" );


	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--receiverheight_km", "Receiver height in km [0.0]" );

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km", 1000.0 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxrange_km", 0.01 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--maxrange_km", "Maximum range in km to use for modeling [1000.0]" );

	ps->addParameter( new NCPA::IntegerParameter( "Nrng_steps", 1000 ) );
	ps->addTest( new NCPA::IntegerGreaterThanTest( "Nrng_steps", 0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nrng_steps", "Number of range steps to use [1000]" );

	ps->addParameter( new NCPA::StringParameter( "ground_impedence_model", "rigid" ) );
	test = ps->addTest( new NCPA::StringSetTest( "ground_impedence_model" ) );
	test->addStringParameter( "rigid" );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence_model", "Impedence model to use.  Currently only \"rigid\" is supported. [rigid]" );

	ps->addParameter( new NCPA::StringParameter( "use_attn_file", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--use_attn_file", "File name containing attenuation, to override default Sutherland/Bass [n/a]. Columns are #n# Height(km) Attenuation(np/m)" );

	ps->addParameter( new NCPA::StringParameter( "dispersion_file", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--dispersion_file", "Filename to output the dispersion information" );

	ps->addParameter( new NCPA::FlagParameter( "append_dispersion_file" ) );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Optional Parameters [default]", "--append_dispersion_file", 
		"Append results to dispersion file rather than overwriting" );
	ps->resetParameterIndent();
	ps->addParameter( new NCPA::StringParameter( "filetag", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--filetag", "Tag prepended to output filenames [n/a]" );


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
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth", "singleprop" ) );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth", "Azimuth of propagation, in degrees CW from North [0,360)" );
	ps->resetParameterIndent();

	
	// for multiple propagation, must specify azimuth start, end, and step
	ps->addParameterDescription( "Modes of Operation", "--multiprop", "Multiple azimuth propagation.  Requires --azimuth_start, --azimuth_end, and --azimuth_step" );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_start" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_start", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_start", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_start", "multiprop" ) );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_end" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_end", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_end", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_end", "multiprop" ) );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_step" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_step", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_step", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_step", "multiprop" ) );
	
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_start", "Starting azimuth, in degrees CW from North [0,360)" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_end", "Ending azimuth, in degrees CW from North [0,360)" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_step", "Azimuth step, in degrees CW from North [0,360)" );
	ps->resetParameterIndent();


	// Setup flags
	//ps->addUsageLine( "Flags:" );
	ps->addParameter( new NCPA::FlagParameter( "write_2d_tloss" ) );
	ps->addParameterDescription( "Flags", "--write_2d_tloss", "Output 2-D transmission loss to tloss2D.nm" );

	ps->addParameter( new NCPA::FlagParameter( "write_lossless" ) );
	ps->addParameterDescription( "Flags", "--write_lossless", "Output lossless as well as lossy results (i.e. omitting atmospheric attenuation).");

	ps->addParameter( new NCPA::FlagParameter( "write_phase_speeds" ) );
	ps->addParameterDescription( "Flags", "--write_phase_speeds", "Output phase speeds to phasespeeds.nm" );

	ps->addParameter( new NCPA::FlagParameter( "write_modes" ) );
	ps->addParameterDescription( "Flags", "--write_modes", "Output modes to mode_###.nm.  Also implies --write_speeds" );

	ps->addParameter( new NCPA::FlagParameter( "write_atm_profile" ) );
	ps->addParameterDescription( "Flags", "--write_atm_profile", "Output atmospheric profile to atm_profile.nm" );

	ps->addParameter( new NCPA::FlagParameter( "Lamb_wave_BC" ) );
	ps->addParameterDescription( "Flags", "--Lamb_wave_BC", "Use admittance = -1/2*dln(rho)/dz" );

	ps->addParameter( new NCPA::FlagParameter( "turnoff_WKB" ) );
	ps->addParameterDescription( "Flags", "--turnoff_WKB", "Turn off the WKB least phase speed estimation" );

	ps->addParameter( new NCPA::FlagParameter( "wvnum_filter" ) );
	//ps->addUsageLine( "  --wvnum_filter          Use wavenumber filter by phase speed.  Requires:" );
	ps->addParameterDescription( "Flags", "--wvnum_filter", "Use wavenumber filter by phase speed.  Requires --c_min and --c_max" );

	ps->addParameter( new NCPA::FloatParameter( "c_min" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "c_min", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "c_min", "wvnum_filter" ) );
	//ps->addUsageLine( "    --c_min               Minimum phase speed to keep" );
	ps->addParameter( new NCPA::FloatParameter( "c_max" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "c_max", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "c_max", "wvnum_filter" ) );
	//ps->addUsageLine( "    --c_max               Maximum phase speed to keep" );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Flags", "--c_min", "Minimum phase speed to keep" );
	ps->addParameterDescription( "Flags", "--c_max", "Maximum phase speed to keep" );
	ps->resetParameterIndent();

	// Dummy parameters
	ps->addParameter( new NCPA::FlagParameter( "broadband" ) );
	ps->addParameter( new NCPA::FloatParameter( "f_min", 0.0 ) );
	ps->addParameter( new NCPA::FloatParameter( "f_max", 0.0 ) );
	ps->addParameter( new NCPA::FloatParameter( "f_step", 0.0 ) );

	// Footer with file formats and sample commands
	ps->addBlankFooterLine();
	ps->addFooterText("OUTPUT Files:  Format description (column order):");
	ps->addFooterTextVerbatim("  tloss_1d.wnm:                r, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)");
	ps->addFooterTextVerbatim("  tloss_1d.lossless.wnm:       r, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)");
	ps->addFooterTextVerbatim("  tloss_2d.wnm:                r, z, 4*PI*Re(P), 4*PI*Im(P)");
	ps->addFooterTextVerbatim("  tloss_2d.lossless.wnm:       r, z, 4*PI*Re(P), 4*PI*Im(P)");
	ps->addFooterTextVerbatim("  Nby2D_tloss_1d.wnm:          r, theta, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)");
	ps->addFooterTextVerbatim("  Nby2D_tloss_1d.lossless.wnm: r, theta, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)");
	ps->addFooterTextVerbatim("  phasespeeds.wnm:             Mode#, phase speed [m/s], imag(k)");
	ps->addFooterTextVerbatim("  mode_<mode_count>.wnm        z, (Mode amplitude)");
	ps->addFooterTextVerbatim("  dispersion_<freq>.wnm        Contains one line with entries:");
	ps->addFooterTextVerbatim("                               freq, (# of modes), rho(z_src), rho(z_rcv)");
	ps->addFooterTextVerbatim("                               followed for each mode 'i' by quadruples:");
	ps->addFooterTextVerbatim("                               real(k(i)), imag(k(i)), Mode(i)(z_src), Mode(i)(z_rcv)");
	ps->addFooterTextVerbatim("  atm_profile.wnm              z,u,v,t,d,p,c,c_eff");
	ps->addBlankFooterLine();
	ps->addFooterText("Examples (run from 'samples' directory):");
	ps->setFooterIndent( 4 );
	ps->setFooterHangingIndent( 4 );
	ps->setCommandMode( true );
	ps->addFooterText("../bin/WMod --singleprop --atmosfile profile_noheader.dat --atmosheaderfile sampleheader.dat --azimuth 90 --freq 0.1");
	ps->addBlankFooterLine();
	ps->addFooterText("../bin/WMod --singleprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat --azimuth 90 --freq 0.1 --write_2d_tloss --sourceheight_km 60 --receiverheight_km 60");
	ps->addBlankFooterLine();
	ps->addFooterText("../bin/WMod --multiprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat --freq 0.1 --azimuth_start 0 --azimuth_end 360 --azimuth_step 1");
	ps->setFooterHangingIndent( 0 );
	ps->setCommandMode( false );
	ps->resetFooterIndent();

}