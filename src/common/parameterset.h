/*
parameterset.h
@version 1.0

The NCPA Parameter Set library.
This library contains objects and methods for parsing, storing, retrieving, and
validating command-line and file-based parameters/options.  Space- and character-
delimited files are supported, as are single-character and multi-character flags
and options on the command line.  All expected options should be provided to the
ParameterSet before parsing.

Option files should be in the format:

# this is a comment
decision : yes
azimuth : 270.0
test_int = 10
verbose
# This is another comment
		
The parameter key and value are separated by one or more delimiter characters.
Delimiters default to ":= " but can be set as needed.  In particular, removing
the space character from the set of delimiters will allow keys or values to include
spaces.  Characters indicating comments can also be set.

Command-line options can be processed as:
programname -f --verbose --azimuth 270.0 --decision=yes

Use of a double hyphen indicates a multi-character key.  If the type of the parameter
expects a value, the next word in the command line is interpreted as the value unless
a delimiter is present, in which case the word is split at the delimiter.  Use of a
single hyphen indicates a single-character flag; if multiple characters are provided, 
such as in '-verbose' as opposed to '--verbose', each letter is treated as an individual
flag.  No values are parsed in these cases.  Any words in the command line that are not
interpreted as options or associated values are preserved and may be retrieved after
parsing.  The original argc and argv values are not modified by this process.

By default the parser runs in strict mode, and will throw an invalid_argument exception
if it encounters a parameter that has not already been defined.  If strict mode is
disabled, unrecognized options will be interpreted as string parameters (if a value is
provided) or flags (if no value is provided or, on the command line, if the next word
can be interpreted as an option or flag, i.e. it begins with - or --).

Validation tests can be added to test parsed values.  Tests are designed to 'fail open',
that is they will default to passing if the value is not set.  There exists a test
specifically for whether a parameter is set or not.  Multiple tests can be designated for
each value, so for example a value indicating an angle can be required, required to be
a positive number, and required to be less than or equal to 360.




Example code:
#include "parameterset.h"

using namespace std;
using namespace NCPA;

int main( int argc, char **argv ) {

	// initialize the master ParameterSet object
	ParameterSet *ps = new ParameterSet();

	// set the delimiter characters, comment characters, and disable strict mode
	ps->setDelimiters( "=|:" );  // default is ":= "
	ps->setComments( "#%!" );    // default is "#"
	ps->setStrict( false );

	// Setup involves telling the ParameterSet what options and flags to expect,
	// what usage text to display, and/or what validation tests should be run on
	// the parameter after it is read

	// This parameter is required, and takes an integer argument that must be positive
	ps->addParameter( new IntegerParameter( "testint" ) );
	ps->addTest( "testint", PARAMETER_TEST_INTEGER_POSITIVE );
	ps->addTest( "testint", PARAMETER_TEST_REQUIRED );

	// This parameter is optional, has a default value of 180.0, and must be in the range
	// [0.0,360.0).
	ps->addParameter( new FloatParameter( "azimuth", 180.0 ) );
	ParameterTest *test = ps->addTest( "azimuth", PARAMETER_TEST_FLOAT_GREATER_THAN_OR_EQUAL );
	test->addFloatParameter( 0.0 );
	test = ps->addTest( "azimuth", PARAMETER_TEST_FLOAT_LESS_THAN );
	test->addFloatParameter( 360.0 );

	// This flag is optional
	ps->addParameter( new FlagParameter( "verbose" ) );

	// This parameter takes a string argument that is optional, but must be one of a set 
	// of values if it is present
	ps->addParameter( new StringParameter( "decision" ) );
	test = ps->addTest( "decision", PARAMETER_TEST_STRING_SET );
	test->addStringParameter( "yes" );
	test->addStringParameter( "no" );
	test->addStringParameter( "maybe" );

	// add help text for the user
	ps->addUsageLine( "  --testint <val>  Test integer value (required)" );
	ps->addUsageLine( "  --azimuth <val>  Floating point value in [0.0,360)" );
	ps->addUsageLine( "  --verbose        Turns verbose mode on" );

	try {
		ncomm = ps->parseFile( "test.options" );
		if (ncomm == 0) {
			cout << "No file options read" << endl;
		}
		ncomm = ps->parseCommandLine( argc, argv );
		if (ncomm == 0) {
			cout << "No command line options read" << endl;
		}
	} catch (invalid_argument &err) {
		cout << "Error parsing arguments: " << err.what() << endl;
		exit( 1 );
	}

	// print a summary of the parameters
	cout << "Parameters:" << endl << endl;
	ps->printParameters( cout );

	// run all of the validation tests
	cout << endl << "Validation:" << endl;
	if (ps->validate()) {
		cout << "All tests passed!" << endl;
	} else {
		ps->printFailedTests( cout );
		ps->printUsage( cout );
	}

	// check if the user wants verbose mode
	if ( ps->getParameter( "verbose" )->wasFound() ) {
		cout << "Verbose mode!" << endl;
	}

	// See if the azimuth value was supplied or default
	if ( ps->getParameter( "azimuth" )->wasFound() ) {
		cout << "Azimuth was user-supplied as " 
			 << ps->getParameter( "azimuth" )->getFloatValue() << endl;
	} else {
		cout << "Azimuth has default value of "
			 << ps->getParameter( "azimuth" )->getStringValue() << endl;
	}
*/

#ifndef _NCPA_PARAMETERSET_H_
#define _NCPA_PARAMETERSET_H_


// user-controlled formatting
#ifndef DEFAULT_TEXT_WIDTH
#define DEFAULT_TEXT_WIDTH 80
#endif

#ifndef DEFAULT_PARAMETER_FIRST_COLUMN_WIDTH
#define DEFAULT_PARAMETER_FIRST_COLUMN_WIDTH 26
#endif

#ifndef DEFAULT_HEADER_INDENT
#define DEFAULT_HEADER_INDENT 0
#endif

#ifndef DEFAULT_FOOTER_INDENT
#define DEFAULT_FOOTER_INDENT 0
#endif

#ifndef DEFAULT_PARAMETER_INDENT
#define DEFAULT_PARAMETER_INDENT 2
#endif







#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>

namespace NCPA {
	
	/**
	  * These indicate tests that do not depend on the type of the tested
	  * value.
	  */
	enum PARAMETER_TEST_TYPE : unsigned int {
	
		// This option/flag must be present.
		PARAMETER_TEST_REQUIRED,
	
		
		// This option/flag must be present if at least one of a set of 
		// other options is present. 	
		PARAMETER_TEST_REQUIRED_IF,
	 
		// Designates a group of options, one and only one of which
		// must be present.
		PARAMETER_TEST_RADIO_BUTTON,
	
		// Integer that must be > 0
		PARAMETER_TEST_INTEGER_POSITIVE,
	
		/** Integer that must be < 0 */
		PARAMETER_TEST_INTEGER_NEGATIVE,
	
		/** Integer that must be > another integer */
		PARAMETER_TEST_INTEGER_GREATER_THAN,
	
		/** Integer that must be >= another integer */
		PARAMETER_TEST_INTEGER_GREATER_THAN_OR_EQUAL,
	
		/** Integer that must be < another integer */
		PARAMETER_TEST_INTEGER_LESS_THAN,
	
		/** Integer that must be <= another integer */
		PARAMETER_TEST_INTEGER_LESS_THAN_OR_EQUAL,
	
		/** Integer that must == 0 */
		PARAMETER_TEST_INTEGER_ZERO,
	
		/** Integer that must != 0 */
		PARAMETER_TEST_INTEGER_NONZERO,
	
		/** Integer that must == another integer */
		PARAMETER_TEST_INTEGER_EQUAL,
	
		/** Integer that must != another integer */
		PARAMETER_TEST_INTEGER_NOT_EQUAL,
	
		/** Double that must > 0.0 */
		PARAMETER_TEST_FLOAT_POSITIVE,
	
		/** Double that must < 0.0 */
		PARAMETER_TEST_FLOAT_NEGATIVE,
	
		/** Double that must > another double */
		PARAMETER_TEST_FLOAT_GREATER_THAN,
	
		/** Double that must >= another double */
		PARAMETER_TEST_FLOAT_GREATER_THAN_OR_EQUAL,
	
		/** Double that must < another double */
		PARAMETER_TEST_FLOAT_LESS_THAN,
	
		/** Double that must <= another double */
		PARAMETER_TEST_FLOAT_LESS_THAN_OR_EQUAL,
	
		/** Double that must == another double (standard floating point caveats apply) */
		PARAMETER_TEST_FLOAT_EQUAL,
	
		/** Double that must != another double (standard floating point caveats apply) */
		PARAMETER_TEST_FLOAT_NOT_EQUAL,
	
		/** Double that must == 0.0 (standard floating point caveats apply) */
		PARAMETER_TEST_FLOAT_ZERO,
	
		/** Double that must != 0.0 (standard floating point caveats apply) */
		PARAMETER_TEST_FLOAT_NONZERO,		
	
		/** string .size() must be >= an integer */
		PARAMETER_TEST_STRING_MINIMUM_LENGTH,
	
		/** string .size() must be <= an integer */
		PARAMETER_TEST_STRING_MAXIMUM_LENGTH,
	
		/** string must match one of a set of strings */
		PARAMETER_TEST_STRING_SET
	};
	

	/**
	  * @brief Abstract base class for command-line and/or file-based parameters.
	  */
	class GenericParameter {
	public:
		virtual ~GenericParameter();

		/** @brief true if the parameter requires an argument, false otherwise */
		virtual bool needsArgument() const = 0;

		/** @brief Parse the argument string as appropriate */
		virtual void parseArgument( const std::string& arg ) = 0;

		/** 
		  * @brief Indicate that the parameter was actually found, and not just set to default value.
		  */
		virtual void setFound( bool tf );

		/** @brief Returns whether the parameter has been found or not. */
		virtual bool wasFound() const;

		/** @brief Whether the parameter has a valid value, which may be default. */
		virtual bool isValid() const;

		/** @brief Return the key string of the parameter */
		virtual std::string getKey() const;

		/** @brief true if the parameter is a Null parameter, false otherwise */
		virtual bool isNull() const;

		/** @brief Get the value of the parameter as an integer. */
		virtual int getIntegerValue() const = 0;

		/** @brief Get the value of the parameter as a floating-point number. */
		virtual double getFloatValue() const = 0;

		/** @brief Get the value of the parameter as a string. */
		virtual std::string getStringValue() const = 0;

		/** @brief Get the value of the parameter as a bool. */
		virtual bool getBoolValue() const = 0;

		/** @brief Print a brief description of the parameter. */
		virtual std::string description() const = 0;

		/** @brief Print the current status of the parameter, including any associated tests. */
		virtual std::string status() const = 0;

	protected:
		std::string _key;
		bool _found;
		bool _valid;
	};

	/**
	  * @brief A type of parameter indicating a null value.
	  * Typically this parameter is returned when requesting a parameter that has 
	  * not been defined.
	  */
	class NullParameter : public GenericParameter {
	public:
		NullParameter( std::string key );
		bool needsArgument() const;
		void parseArgument( const std::string& arg );
		bool isNull() const;
		
		// these will always throw invalid_argument
		int getIntegerValue() const;
		double getFloatValue() const;
		std::string getStringValue() const;
		bool getBoolValue() const;

		// these will return an empty string
		std::string description() const;
		std::string status() const;
	};

	/**
	  * @brief A type of parameter representing an integer value.
	  */
	class IntegerParameter : public GenericParameter {
	public:
		IntegerParameter( std::string key );
		IntegerParameter( std::string key, int defaultValue );
		void setValue( int newval );
		bool needsArgument() const;
		void parseArgument( const std::string& arg );
		int getIntegerValue() const;
		double getFloatValue() const;
		std::string getStringValue() const;
		bool getBoolValue() const;
		std::string description() const;
		std::string status() const;
	protected:
		int _value;
	};
	
	/**
	  * @brief A type of parameter indicating a floating-point value.
	  */
	class FloatParameter : public GenericParameter {
	public:
		FloatParameter( std::string key );
		FloatParameter( std::string key, double defaultValue );
		void setValue( double newval );
		bool needsArgument() const;
		void parseArgument( const std::string& arg );
		int getIntegerValue() const;
		double getFloatValue() const;
		std::string getStringValue() const;
		bool getBoolValue() const;
		std::string description() const;
		std::string status() const;
	protected:
		double _value;
	};

	/**
	@brief A type of parameter indicating a string value.
	*/
	class StringParameter : public GenericParameter {
	public:
		StringParameter( std::string key );
		StringParameter( std::string key, std::string defaultValue );
		void setValue( std::string newval );
		bool needsArgument() const;
		void parseArgument( const std::string& arg );
		int getIntegerValue() const;
		double getFloatValue() const;
		std::string getStringValue() const;
		bool getBoolValue() const;
		std::string description() const;
		std::string status() const;
	protected:
		std::string _value;
	};
	
	/**
	@brief A type of parameter indicating a flag with no associated value.
	*/
	class FlagParameter : public GenericParameter {
	public:
		FlagParameter( std::string key );
		FlagParameter( std::string key, bool defaultValue );
		void setValue( bool newval );
		bool needsArgument() const;
		void parseArgument( const std::string& arg );
		int getIntegerValue() const;
		double getFloatValue() const;
		std::string getStringValue() const;
		bool getBoolValue() const;
		bool isValid() const;
		std::string description() const;
		std::string status() const;
	protected:
		bool _value;
	};

	/**
	Extension of a std::vector< NCPA::GenericParameter * > with a custom search method.
	Should not be directly instantiated.
	@brief Extension of std::vector< NCPA::GenericParameter * >
	*/
	class ParameterVector : public std::vector< NCPA::GenericParameter * > {
	public:
		GenericParameter * findParameter( const std::string& key ) const;
	};

/**
	  * Abstract base class for validation criteria.  See subclass descriptions
	  * for general usage, you won't instantiate this class directly.
	  */
	class ParameterTest {
	public:

		/**
		Default constructor.
		@brief Default constructor.
		*/
		ParameterTest();

		/**
		Initialization constructor.
		@brief Initialization constructor for inherited classes to call.
		*/
		ParameterTest( std::string option_name, bool starts_ready );

		/**
		  * Destructor.  Cleans up any dynamically allocated memory.
		  */
		virtual ~ParameterTest();
	
		/**
		  * Run the validation checks specified for this option
		  * @param opts A pointer to the AnyOption object that has ingested
		  *             the command line and file options
		  * @return true if the test passes, false otherwise
		  */
		virtual bool validate( const ParameterVector& param ) = 0;
	
		/**
		  * A text description of the test that is to be run.
		  * @return A std::string containing a description of the test
		  */
		virtual std::string description() const = 0;
	
		/**
		  * A text description of why the test failed, if applicable.
		  * @return A std::string containing a description of the failure
		  */
		virtual std::string failureMessage() const = 0;
	
		/**
		  * The name of the parameter to be checked
		  * @return A std::string containing the parameter name
		  */
		virtual std::string optionName() const;
	
		/**
		  * Add an integer parameter to the test, if applicable.
		  *
		  * Does nothing unless overridden.
		  * @param param	The parameter to add to the test
		  */
		virtual void addIntegerParameter( int param );
	
		/**
		  * Add a double parameter to the test, if applicable.
		  *
		  * Does nothing unless overridden.
		  * @param param	The parameter to add to the test
		  */
		virtual void addFloatParameter( double param );
	
		/**
		  * Add a string parameter to the test, if applicable.
		  *
		  * Does nothing unless overridden.
		  * @param param	The parameter to add to the test
		  */
		virtual void addStringParameter( const std::string param );
	
		/**
		  * Indicates if the test is ready to be run (i.e. any necessary 
		  * parameters have been supplied).
		  *
		  * @return true if the test can be run meaningfully, false otherwise
		  */
		virtual bool ready() const;
	
		/**
		  * Returns the value the test is checking for, as a string.
		  *
		  * @return the test value in string form, implementation-dependent
		  */
		virtual std::string valueString() const = 0;
	
	protected:
	
		/**
		  * The option name.
		  */
		std::string _optName;
	
		/**
		  * All information has been provided and the test can be run.
		  */
		bool _ready;
	};







	/**
	Handles the specification, ingestion, and tracking of user-supplied flags and input
	arguments, including help text and argument validation.
	@brief Handles input arguments and parameters
	*/
	class ParameterSet {
	public:

		/**
		Default constructor
		@brief Default constructor
		*/
		ParameterSet();

		/**
		Destructor.
		@brief Destructor.
		*/
		~ParameterSet();

		/**
		Changes the delimiters to use for splitting up parameter file entries into keys and
		values.  Multiple consecutive delimiters will be treated as one.
		@brief Set the parameter file key/value delimiters.
		@param newdelims A string containing the new delimiter character set.
		*/
		void setDelimiters( std::string newdelims );

		/**
		Changes the character(s) to be used to indicate comment lines in parameter files.
		@brief Changes the character(s) to be used to indicate comment lines in parameter files.
		@param newcomms A string containing the new comment character set.
		*/
		void setComments( std::string newcomms );

		/**
		Add a new parameter to the expected parameter set.  Usually called using the constructor
		for the specific parameter type in question, e.g.

		ps->addParameter( new NCPA::FlagParameter( "flagname" ) );

		@brief Add a new parameter to the set
		@param newParam A pointer to a specific instance of a GenericParameter subclass
		*/
		void addParameter( GenericParameter *newParam );

		/**
		Turns strict mode on or off.  When in strict mode, the ParameterSet will thrown an exception
		if it parses a parameter that it was not told to expect.  Otherwise, it will silently add the 
		parameter as a FlagParameter.
		@brief Set strict mode on or off
		@param tf The new state for strict mode
		*/
		void setStrict( bool tf );

		// Messages to output to the user
		//void addUsageLine( const std::string& line );
		void setHeaderIndent( unsigned int newindent );
		void resetHeaderIndent();
		void setHeaderHangingIndent( unsigned int newindent );
		void addBlankHeaderLine();
		void addHeaderText( const std::string& text );
		void addHeaderTextVerbatim( const std::string& text );

		void setFooterIndent( unsigned int newindent );
		void resetFooterIndent();
		void setFooterHangingIndent( unsigned int newindent );
		void addBlankFooterLine();
		void addFooterText( const std::string& text );
		void addFooterTextVerbatim( const std::string& text );
		
		void setParameterIndent( unsigned int newindent );
		void resetParameterIndent();
		void addParameterDescription( const std::string& section, const std::string& param, 
			const std::string &description,
			unsigned int firstcolumnwidth = DEFAULT_PARAMETER_FIRST_COLUMN_WIDTH );
		
		void setCommandMode( bool tf );
		void setTextWidth( unsigned int newWidth );
		void printUsage( std::ostream& os = std::cout ) const;

		// Ingest from command line or file
		unsigned int parseCommandLine( unsigned int argc, char **argv );
		unsigned int parseFile( std::string filename );
		std::vector< std::string > getUnparsedOptions() const;

		// Testing of parsed parameters
		ParameterTest * addTest( const std::string &key, PARAMETER_TEST_TYPE option_type );
		ParameterTest * addTest( ParameterTest *newTest );
		bool validate();
		void printFailedTests( std::ostream& os = std::cerr ) const;
		std::vector< ParameterTest * > getTests( std::string key ) const;
		
		void printParameters( bool printTests = false, std::ostream& os = std::cout ) const;

		GenericParameter *getParameter( std::string key );
		int getInteger( std::string key ) const;
		double getFloat( std::string key ) const;
		std::string getString( std::string key ) const;
		bool getBool( std::string key ) const;     // returns bool version of parameter
		bool wasFound( std::string key ) const;    // returns true if the value was specified, not left default
		void removeParameter( std::string key );
		
	protected:
		ParameterVector _params;
		std::string _delims, _comments;
		std::vector< std::string > _unparsed;
		std::vector< std::string > _headerLines, _footerLines, _sections;
		std::map< std::string, std::ostringstream * > _descriptionLines;
		std::vector< NCPA::ParameterTest * > _tests, _failed_tests;
		bool _strict, _commandMode;
		unsigned int headerIndent_, footerIndent_, parameterIndent_, 
			maxWidth_, headerHangingIndent_, footerHangingIndent_;


		bool isLongOption_( std::string ) const;
		bool isShortOption_( std::string ) const;
		unsigned int processShortOption_( int argc, char **argv, unsigned int index );
		unsigned int processLongOption_( int argc, char **argv, unsigned int index );
		void processSingleOption_( std::string key );
		void processDoubleOption_( std::string key, std::string value );

		void formatText_( std::vector< std::string > &holder, const std::string& text, 
			unsigned int indent, unsigned int hanging_indent, unsigned int maxwidth );
		void addSpaces_( std::ostringstream *oss, unsigned int spaces );
		
	};





	/** Test whether an option is present */
	class RequiredTest : public ParameterTest {
	public:
		RequiredTest( const std::string& option_name );
		bool validate( const ParameterVector& param );
		std::string description() const;
		std::string failureMessage() const;
		std::string valueString() const;
	};
	
	/** Test whether an option is present if another option is also set */
	class RequiredIfOtherIsPresentTest : public ParameterTest {
	public:
		RequiredIfOtherIsPresentTest( const std::string& option_name );
		RequiredIfOtherIsPresentTest( const std::string& option_name,
			const std::string& prereq );
		RequiredIfOtherIsPresentTest( const std::string& option_name, 
			unsigned int nPrereqs, std::string *prereqs );
		RequiredIfOtherIsPresentTest( const std::string& option_name,
			const std::vector< std::string > prereq_vector );
		bool validate( const ParameterVector& param );
		std::string description() const;
		std::string failureMessage() const;
		std::string valueString() const;
		void addStringParameter( const std::string param );
		bool ready() const;
	private:
		std::vector< std::string > _prereqs;
	};
	
	/** Test whether one and only one of a set of options or flags is present */
	class RadioButtonTest : public ParameterTest {
	public:
		RadioButtonTest( const std::string& option_name );
		RadioButtonTest( const std::string& option_name, unsigned int nbuttons, 
			std::string *buttons );
		RadioButtonTest( const std::string& option_name, 
			const std::vector< std::string > buttons );
		bool validate( const ParameterVector& param );
		std::string description() const;
		std::string failureMessage() const;
		void addStringParameter( const std::string param );
		std::string valueString() const;
		std::vector< std::string > lastMatched() const;
		bool ready() const;
	private:
		std::vector< std::string > _buttons;
		std::vector< std::string > _matched;
	};

	/** Test whether an integer option is greater than a parameter */
	class IntegerGreaterThanTest : public ParameterTest {
	public:
		IntegerGreaterThanTest( const std::string& option_name );
		IntegerGreaterThanTest( const std::string& option_name, int comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value, _testedValue;
	};

	/** Test whether an integer option is greater than or equal to a parameter */
	class IntegerGreaterThanOrEqualToTest : public ParameterTest {
	public:
		IntegerGreaterThanOrEqualToTest( const std::string& option_name );
		IntegerGreaterThanOrEqualToTest( const std::string& option_name, int comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value, _testedValue;
	};

	/** Test whether an integer option is less than a parameter */
	class IntegerLessThanTest : public ParameterTest {
	public:
		IntegerLessThanTest( const std::string& option_name );
		IntegerLessThanTest( const std::string& option_name, int comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value, _testedValue;
	};
	
	// Test whether an integer option is less than
	class IntegerLessThanOrEqualToTest : public ParameterTest {
	public:
		IntegerLessThanOrEqualToTest( const std::string& option_name );
		IntegerLessThanOrEqualToTest( const std::string& option_name, int comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value, _testedValue;
	};
	
	// Test whether an integer option is equal to
	class IntegerEqualToTest : public ParameterTest {
	public:
		IntegerEqualToTest( const std::string& option_name );
		IntegerEqualToTest( const std::string& option_name, int comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		int _value, _testedValue;
	};
	
	// Test whether an integer option is not equal to
	class IntegerNotEqualToTest : public ParameterTest {
	public:
		IntegerNotEqualToTest( const std::string& option_name );
		IntegerNotEqualToTest( const std::string& option_name, int comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		int _value, _testedValue;
	};
	
	
	// Test whether a floating point option is greater than
	class FloatGreaterThanTest : public ParameterTest {
	public:
		FloatGreaterThanTest( const std::string& option_name );
		FloatGreaterThanTest( const std::string& option_name, double comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value, _testedValue;
	};
	
	// Test whether a floating point option is greater than or equal to
	class FloatGreaterThanOrEqualToTest : public ParameterTest {
	public:
		FloatGreaterThanOrEqualToTest( const std::string& option_name );
		FloatGreaterThanOrEqualToTest( const std::string& option_name, double comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value, _testedValue;
	};

	// Test whether a floating point option is less than
	class FloatLessThanTest : public ParameterTest {
	public:
		FloatLessThanTest( const std::string& option_name );
		FloatLessThanTest( const std::string& option_name, double comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value, _testedValue;
	};
	
	// Test whether a floating point option is less than or equal to
	class FloatLessThanOrEqualToTest : public ParameterTest {
	public:
		FloatLessThanOrEqualToTest( const std::string& option_name );
		FloatLessThanOrEqualToTest( const std::string& option_name, double comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value, _testedValue;
	};
	
	// Test whether a floating point option is equal to (standard caveats apply)
	class FloatEqualToTest : public ParameterTest {
	public:
		FloatEqualToTest( const std::string& option_name );
		FloatEqualToTest( const std::string& option_name, double comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value, _testedValue;
	};
	
	// Test whether a floating point option is not equal to (standard caveats apply)
	class FloatNotEqualToTest : public ParameterTest {
	public:
		FloatNotEqualToTest( const std::string& option_name );
		FloatNotEqualToTest( const std::string& option_name, double comparison );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value, _testedValue;
	};
	
	// Test whether the length of a string is at least N characters
	class StringMinimumLengthTest : public ParameterTest {
	public:
		StringMinimumLengthTest( const std::string& option_name );
		StringMinimumLengthTest( const std::string& option_name, size_t minlength );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		size_t _value;
		std::string _testedValue;
	};
	
	// Test whether the length of a string is at most N characters
	class StringMaximumLengthTest : public ParameterTest {
	public:
		StringMaximumLengthTest( const std::string& option_name );
		StringMaximumLengthTest( const std::string& option_name, size_t maxlength );
		bool validate( const ParameterVector& param ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		size_t _value;
		std::string _testedValue;
	};
	
	
	// test whether a string is a member of a set or not
	class StringSetTest : public ParameterTest {
	public:
		StringSetTest( const std::string& option_name );
		StringSetTest( const std::string& option_name, unsigned int nSet, 
			std::string *choices );
		StringSetTest( const std::string& option_name, 
			std::vector< std::string > choices );
		bool validate( const ParameterVector& param );
		std::string description() const;
		std::string failureMessage() const;
		void addStringParameter( const std::string choice_name );
		std::string valueString() const;
		bool ready() const;
	private:
		std::vector< std::string > _choices;
		std::string _testedValue;
	};


}




#endif