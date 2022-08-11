#include "ncpaprop_common.h"
#include "epade_pe_parameters.h"
#include <string>
#include <iostream>
#include <vector>


void NCPA::configure_epade_solver( NCPA::EPadeSolver *solver,
	NCPA::ParameterSet *param ) {

	NCPA::units_t u_km = NCPA::Units::fromString("km"),
				  u_m  = NCPA::Units::fromString("m"),
				  u_K  = NCPA::Units::fromString("K");
	double tempf1, tempf2, tempf3;
	int    tempi;
	std::string temps;
	std::vector<double> tempv;

	// simple values - use defaults if unspecified
	solver->setMaximumRange( param->getFloat( "maxrange_km" ), u_km );
	if (param->wasFound("maxheight_km")) {
		solver->setMaximumHeight( param->getFloat( "maxheight_km" ), u_km );
	}
	if (param->wasFound("sourceheight_km")) {
		solver->setSourceHeight( param->getFloat("sourceheight_km", u_km ) );
	}
	if (param->wasFound("receiverheight_km")) {
		solver->setReceiverHeight( param->getFloat("receiverheight_km", u_km ) );
	}

	solver->setHeightResolution( param->getFloat("dz_m"), u_m );
	solver->setRangeSteps( param->getInteger("Nrng_steps") );
	solver->setPadeOrder( param->getInteger( "npade" ) );
	solver->setFileTag( param->getString( "filetag" ) );

	// only do these if they're specifically there
	if (param->wasFound( "ground_impedence_imag") ||
		param->wasFound( "ground_impedence_real") ) {

		solver->setUserGroundImpedence(
			param->getFloat( "ground_impedence_real" ),
			param->getFloat( "ground_impedence_imag" ) );
	}
	if (param->wasFound( "groundheight_km" ) ) {
		solver->setUserGroundHeight( param->getFloat( "groundheight_km"),
			u_km );
	}

	// atmosphere
	if (param->wasFound( "atmosfile" )) {
		solver->setAtmosphere(
			NCPA::AtmosphereType::STRATIFIED_2D,
			param->getString( "atmosfile" ),
			param->getString( "atmosheaderfile" ) );
	} else if (param->wasFound( "atmosfile2d" ) ) {
		solver->setAtmosphere(
			NCPA::AtmosphereType::RANGE_DEPENDENT_2D,
			param->getString( "atmosfile2d" ),
			param->getString( "atmosheaderfile" ) );
	}

	// starter
	temps = param->getString( "starter" );
	if (temps == "self") {
		solver->setStarter( NCPA::StarterType::SELF );
	} else if (temps == "gaussian") {
		solver->setStarter( NCPA::StarterType::GAUSSIAN );
	} else if (temps == "user") {
		solver->setStarter( NCPA::StarterType::USER );
		solver->setStarterFileName( param->getString( "starterfile" ) );
	}

	// azimuth
	if (param->wasFound("singleprop")) {
		solver->setAzimuth( param->getFloat( "azimuth" ) );
	} else if (param->wasFound("multiprop")) {
		double min_az, max_az, step_az;
		min_az = param->getFloat("azimuth_start");
		max_az = param->getFloat("azimuth_end");
		step_az = param->getFloat("azimuth_step");
		while (min_az > max_az) {
			min_az -= 360.0;
		}
		tempv.clear();
		while(min_az <= max_az) {
			tempv.push_back( min_az );
			min_az += step_az;
		}
		solver->setAzimuth( tempv );
		tempv.clear();
	}

	// topography
	if (param->wasFound("topo")) {
		solver->setUseTopography( true );
		solver->setTopographyFileName( param->getString( "topofile" ) );
	}

	// attenuation
	if (param->wasFound("attnfile")) {
		solver->setUserAttenuationUsed( true );
		solver->setAttenuationFileName( param->getString( "attnfile" ) );
	}

	// frequency
	solver->setFrequency( param->getFloat( "freq" ) );

	// optional turbulence
	if (param->wasFound("turbulence")) {
		solver->setTurbulenceScale(
			param->getFloat( "turbulence_scale_m"), u_m );
		solver->setTurbulenceSpectrumSize( param->getInteger( "n_turbulence" ) );
		solver->setTurbulenceTemperatureFactor(
			param->getFloat( "turbulence_t_factor" ) );
		solver->setTurbulenceVelocityFactor(
			param->getFloat( "turbulence_v_factor" ) );
		solver->setTurbulenceFileName(
			param->getString( "turbulence_file" ) );
		solver->setUseTurbulence( true );
	}

	// flags
	solver->setWrite2DTransmissionLoss( param->wasFound( "write_2d_tloss" ) );
	solver->setLossless( param->wasFound( "lossless" ) );
	solver->setUseTopLayer( !param->wasFound( "disable_top_layer" ) );
	solver->setWriteStarter( param->wasFound( "write_starter" ) );
	solver->setWriteTopography( param->wasFound( "write_topography" ) );
	solver->setWriteAtmosphere( param->wasFound( "write_atm_profile" ) );
}


void NCPA::configure_epade_pe_parameter_set( NCPA::ParameterSet *ps ) {

	//NCPA::ParameterTest *test = NULL;

	// general configuration
	// set up expected commands
	ps->setStrict( true );
	ps->setComments( "#%" );
	ps->setDelimiters( "=: " );

	// Add header instructions
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addHeaderTextVerbatim("|                             NCPA Infrasound                              |");
	ps->addHeaderTextVerbatim("|                         ePade Parabolic Equation                         |");
	ps->addHeaderTextVerbatim("|         Single Frequency - Effective Sound Speed Approximation           |");
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addBlankHeaderLine();
	ps->addHeaderText("By default the program computes both the 1D (at the ground) and 2D transfer function and saves the data to 2 files:" );
	ps->setHeaderIndent( 4 );
	ps->addHeaderText("file tloss_1d.pe - transfer function at the ground" );
	ps->addHeaderText("file tloss_2d.pe - full 2-D transfer function field" );
	ps->resetHeaderIndent();
	ps->addBlankHeaderLine();
	ps->addHeaderText("The options below can be specified in a colon-separated file \"epade_pe.param\" or at the command line. Command-line options override file options.");

	// Parameter descriptions
	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameter( new FlagParameter( "h" ) );
	ps->addParameterDescription( "Options Control", "--help", "Prints help test" );

	ps->addParameter( new StringParameter( "paramfile", "epade_pe.param") );
	ps->addParameterDescription( "Options Control", "--paramfile", "Parameter file name [epade_pe.param]" );

	ps->addParameter( new FlagParameter( "printparams" ) );
	ps->addParameterDescription( "Options Control", "--printparams", "Print parameter summary to screen" );

	// Atmosphere
	std::string atmosphere_types[ 2 ] = { "atmosfile", "atmosfile2d" };
	ps->addParameter( new NCPA::StringParameter( atmosphere_types[ 0 ] ) );
	ps->addParameter( new NCPA::StringParameter( atmosphere_types[ 1 ] ) );
	ps->addTest( new NCPA::RadioButtonTest( "atmosphere_type", 2, atmosphere_types ) );
	ps->addParameterDescription( "Atmosphere", "--atmosfile", "1-D atmospheric profile filename" );
	ps->addParameterDescription( "Atmosphere", "--atmosfile2d", "2-D atmospheric summary filename (see manual)" );

	// Required parameters
	ps->addParameter( new NCPA::FloatParameter( "freq" ) );
	ps->addTest( new NCPA::RequiredTest( "freq" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "freq", 0.0 ) );
	ps->addParameterDescription( "Required Parameters", "--freq", "Frequency of analysis (Hz)" );

	ps->addParameter( new NCPA::StringParameter( "starter", "self" ) );
	ps->addTest( new NCPA::RequiredTest( "starter" ) );
	NCPA::ParameterTest *test = ps->addTest( new NCPA::StringSetTest( "starter" ) );
	test->addStringParameter( "self" );
	test->addStringParameter( "gaussian" );
	test->addStringParameter( "user" );
	ps->addParameterDescription( "Required Parameters", "--starter", "Starter type: one of { self, gaussian, user }" );
	

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km" ) );
	ps->addTest( new NCPA::RequiredTest( "maxrange_km" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxrange_km", 0.01 ) );
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
	ps->addParameter( new NCPA::StringParameter( "atmosheaderfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--atmosheaderfile", "External header file, overrides internal header [None]" );
	
	ps->addParameter( new NCPA::IntegerParameter( "npade", 4 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "npade", 3 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--npade", "Number of Pade coefficients to use [4]" );

	ps->addParameter( new NCPA::FloatParameter( "dz_m", -1.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--dz_m", "Altitude resolution in meters [automatic]" );

	ps->addParameter( new NCPA::FloatParameter( "maxheight_km", 150.0 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxheight_km", 0.1 ));
	ps->addParameterDescription( "Optional Parameters [default]", "--maxheight_km", "Maximum height of analysis in km [150.0, or max height of atmosphere]" );

	ps->addParameter( new NCPA::FloatParameter( "sourceheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--sourceheight_km", "Source height in km [ground]" );

	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--receiverheight_km", "Receiver height in km [ground]" );

	ps->addParameter( new NCPA::FloatParameter( "groundheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--groundheight_km", "Ground height in km [Z0 parameter in profile, or 0.0]" );

	ps->addParameter( new NCPA::IntegerParameter( "Nrng_steps", 0 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "Nrng_steps", 0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nrng_steps", "Number of range steps to use [decide internally]" );

	ps->addParameter( new NCPA::FloatParameter( "ground_impedence_real", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence_real", "Real part of ground impedence [rigid ground]" );
	ps->addParameter( new NCPA::FloatParameter( "ground_impedence_imag", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence_imag", "Imaginary part of ground impedence [rigid ground]" );
	ps->addParameter( new NCPA::IntegerParameter( "n_turbulence", 20 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--n_turbulence", "Number of random turbulence phases to compute [20]" );
	ps->addParameter( new NCPA::StringParameter( "turbulence_file" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_file", "File containing 2*n_turbulence numbers in [0,1) range [randomize]" );
	ps->addParameter( new NCPA::FloatParameter( "turbulence_scale_m", 100.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_scale_m", "Turbulence scale (m) [100]" );
	ps->addParameter( new NCPA::FloatParameter( "turbulence_t_factor", 1.0e-10 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_t_factor", "Temperature factor for turbulence spectrum (m^-(2/3)) [1.0e-10]" );
	ps->addParameter( new NCPA::FloatParameter( "turbulence_v_factor", 1.0e-8 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--turbulence_v_factor", "Wind velocity factor for turbulence spectrum (m^-(2/3)) [1.0e-8]" );
	ps->addParameter( new NCPA::StringParameter( "topofile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--topofile", "File name containing topography [n/a]. Columns are #n# Range(km) Elevation(m) by default. #n# Units can be specified using header: #n#  #% r, km #n#  #% z, km" );
	ps->addParameter( new NCPA::StringParameter( "starterfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--starterfile", "File name containing starter [n/a]. Columns are #n# Height(km) RealPart ImaginaryPart" );
	ps->addParameter( new NCPA::StringParameter( "attnfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--attnfile", "File name containing attenuation, to override default Sutherland/Bass [n/a]. Columns are #n# Height(km) Attenuation(np/m)" );
	ps->addParameter( new NCPA::StringParameter( "filetag", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--filetag", "Tag prepended to output filenames [n/a]" );

	// Setup flags
	//ps->addUsageLine( "Flags:" );
	ps->addParameter( new NCPA::FlagParameter( "write_2d_tloss" ) );
	ps->addParameterDescription( "Flags", "--write_2d_tloss", "Output 2-D transfer function to tloss_2d.pe" );
	ps->addParameter( new NCPA::FlagParameter( "write_starter" ) );
	ps->addParameterDescription( "Flags", "--write_starter", "Output starter to starter.pe" );
	ps->addParameter( new NCPA::FlagParameter( "write_topography" ) );
	ps->addParameterDescription( "Flags", "--write_topography", "Output interpolated topography to topography.pe" );
	ps->addParameter( new NCPA::FlagParameter( "write_atm_profile" ) );
	ps->addParameterDescription( "Flags", "--write_atm_profile", "Output atmospheric profile at the source to atm_profile.nm" );
	ps->addParameter( new NCPA::FlagParameter( "lossless" ) );
	ps->addParameterDescription( "Flags", "--lossless", "Ignore atmospheric attenuation" );
	ps->addParameter( new NCPA::FlagParameter( "topo" ) );
	ps->addParameterDescription( "Flags", "--topo", "Use topography.  Requires presence of 'Z0' parameter in atmospheric files" );
	ps->addParameter( new NCPA::FlagParameter( "turbulence" ) );
	ps->addParameterDescription( "Flags", "--turbulence", "Include turbulence." );

	// hidden parameters
	ps->addParameter( new NCPA::FlagParameter( "disable_top_layer" ) );
	ps->addParameter( new NCPA::FlagParameter( "broadband" ) );
	//ps->addParameterDescription( "Flags", "--broadband", "Calculate at multiple frequencies" );


	// Footer with file formats and sample commands
	ps->addBlankFooterLine();
	ps->addFooterText("OUTPUT Files:  Format description (column order):");
	ps->addFooterTextVerbatim("  tloss_1d.pe:                 r (km), az (deg), TF (real), TF (imag)");
	ps->addFooterTextVerbatim("  tloss_multiplot.pe:          r (km), az (deg), TF (real), TF (imag)");
	ps->addFooterTextVerbatim("  tloss_2d.pe:                 r, z, TF (real), TF (imag)");
	ps->addFooterTextVerbatim("  atm_profile.pe:              z,u,v,w,t,d,p,c,c_eff");
	ps->addFooterTextVerbatim("  starter.pe:                  z, starter (real), starter(imag)" );
	ps->addFooterTextVerbatim("  topography.pe:               az, r, z0" );
	ps->addBlankFooterLine();
	ps->addFooterText("Examples (run from 'samples' directory):");
	ps->setFooterIndent( 4 );
	ps->setFooterHangingIndent( 4 );
	ps->setCommandMode( true );
	ps->addBlankFooterLine();

	// example 1
	ps->addFooterText("../bin/ePape --singleprop --starter self --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addBlankFooterLine();

	// example 2
	ps->addFooterText("../bin/ePape --singleprop --starter self --atmosfile2d toy_profile_2d_summary.dat --freq 0.5 --azimuth 90 --maxrange_km 1800 --write_2d_tloss" );
	ps->addBlankFooterLine();

	// example 3
	ps->addFooterText("../bin/ePape --multiprop --starter self --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 --maxrange_km 1000");
	ps->addBlankFooterLine();

	// example 4
	ps->addFooterText("../bin/ePape --singleprop --topo --starter self --atmosfile2d toy_profile_2d_summary_topo.dat --freq 0.5 --azimuth 90 --maxrange_km 1200 --write_2d_tloss --write_topography");
	ps->setFooterHangingIndent( 0 );
	ps->setCommandMode( false );
	ps->resetFooterIndent();
}

