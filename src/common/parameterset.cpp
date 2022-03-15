#include "parameterset.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <cctype>
#include <sstream>
#include "util.h"


/*
*********************************************************************************
Code for ParameterSet class
*********************************************************************************
*/
NCPA::ParameterSet::ParameterSet() 
	: _delims{ ":= " }, _comments{ "#" }, _strict{ true }, _commandMode{ false },
	  headerIndent_{ DEFAULT_HEADER_INDENT }, footerIndent_{ DEFAULT_FOOTER_INDENT },
	  parameterIndent_{ DEFAULT_PARAMETER_INDENT }, maxWidth_{ DEFAULT_TEXT_WIDTH },
	  headerHangingIndent_{ 0 }, footerHangingIndent_{ 0 } {}


NCPA::ParameterSet::~ParameterSet() {
	for ( std::vector< NCPA::GenericParameter * >::iterator it = _params.begin(); 
			it != _params.end(); ++it ) {
		delete *it;
	}
	_params.clear();
	for (std::vector< NCPA::ParameterTest * >::iterator it2 = _tests.begin();
			it2 != _tests.end(); ++it2 ) {
		delete *it2;
	}

	_tests.clear();
	_failed_tests.clear();
	_unparsed.clear();
	//_usage.clear();
	_headerLines.clear();
	_footerLines.clear();

	for (std::vector< std::string >::iterator it3 = _sections.begin(); it3 != _sections.end(); ++it3) {
		delete _descriptionLines.at( *it3 );
	}

	_descriptionLines.clear();
	_sections.clear();
}

void NCPA::ParameterSet::setCommandMode( bool tf ) {

	if (tf != _commandMode) {
		// we're changing, so we need to adjust window width
		if (tf) {
			// entering command mode, window width needs to shrink by 2
			maxWidth_ -= 2;
		} else {
			// leaving command mode, expand window by 2
			maxWidth_ += 2;
		}
	}

	_commandMode = tf;
}

void NCPA::ParameterSet::setTextWidth( unsigned int newWidth ) {
	if (_commandMode) {
		newWidth -= 2;
	}
	maxWidth_ = newWidth;
}

void NCPA::ParameterSet::setHeaderIndent( unsigned int newindent ) {
	headerIndent_ = newindent;
}

void NCPA::ParameterSet::resetHeaderIndent() {
	headerIndent_ = DEFAULT_HEADER_INDENT;
}

void NCPA::ParameterSet::resetFooterIndent() {
	footerIndent_ = DEFAULT_FOOTER_INDENT;
}

void NCPA::ParameterSet::setFooterIndent( unsigned int newindent ) {
	footerIndent_ = newindent;
}

void NCPA::ParameterSet::setParameterIndent( unsigned int newindent ) {
	parameterIndent_ = newindent;
}

void NCPA::ParameterSet::resetParameterIndent() {
	parameterIndent_ = DEFAULT_PARAMETER_INDENT;
}

void NCPA::ParameterSet::setHeaderHangingIndent( unsigned int newindent ) {
	headerHangingIndent_ = newindent;
}

void NCPA::ParameterSet::setFooterHangingIndent( unsigned int newindent ) {
	footerHangingIndent_ = newindent;
}

void NCPA::ParameterSet::addHeaderText( const std::string& text ) {
	formatText_( _headerLines, text, headerIndent_, headerHangingIndent_, maxWidth_ );
}

void NCPA::ParameterSet::addHeaderTextVerbatim( const std::string& text ) {
	std::string tmpStr( text );
	_headerLines.push_back( tmpStr );
}

void NCPA::ParameterSet::addBlankHeaderLine() {
	_headerLines.push_back("");
}

void NCPA::ParameterSet::addBlankFooterLine() {
	_footerLines.push_back("");
}

void NCPA::ParameterSet::addFooterText( const std::string& text ) {
	formatText_( _footerLines, text, footerIndent_, footerHangingIndent_, maxWidth_ );
}

void NCPA::ParameterSet::addFooterTextVerbatim( const std::string& text ) {
	std::string tmpStr( text );
	_footerLines.push_back( tmpStr );
}

void NCPA::ParameterSet::formatText_( std::vector< std::string > &holder, const std::string& text, 
	unsigned int indent, unsigned int hanging_indent, unsigned int maxwidth ) {

	std::string tmpStr;
	std::ostringstream oss;
	//unsigned int indent_i;
	unsigned int hang = 0;

	std::vector< std::string > words = NCPA::split( text );
	addSpaces_( &oss, indent + hang );

	for (std::vector< std::string >::const_iterator it = words.cbegin(); 
		it != words.end(); ++it) {

		// check for specific newline instruction
		if (*it == "#n#") {
			tmpStr = oss.str();
			if (_commandMode) { 
				tmpStr += " \\"; 
			}
			holder.push_back( tmpStr );
			hang = hanging_indent;
			oss.str( "" );
			addSpaces_( &oss, indent + hang );
		} else {

			if ( (oss.str().size() + (*it).size() + 1) > maxwidth ) {
				// adding this word would extend the line too far, so flush it
				tmpStr = oss.str();
				if (_commandMode) { 
					tmpStr += " \\"; 
				}
				holder.push_back( tmpStr );
				hang = hanging_indent;
				oss.str( "" );
				addSpaces_( &oss, indent + hang );
			}

			if (oss.str().size() > indent) {
				oss << " ";
			}
			oss << *it;
		}
	}
	tmpStr = oss.str();
	holder.push_back( tmpStr );
}

void NCPA::ParameterSet::addParameterDescription( const std::string& section, const std::string& param, 
			const std::string &description, unsigned int firstcolumnwidth ) {

	// do we have this header already?
	std::ostringstream *oss, *oss_orig;
/*
	try {
		oss_orig = _descriptionLines.at( section );
	} catch (std::out_of_range &oor) {
		oss_orig = new std::ostringstream("");
		_descriptionLines[ section ] = oss_orig;
	}
*/
	if ( _descriptionLines.count( section ) == 1 ) {
		oss_orig = _descriptionLines[ section ];
	} else {
		oss_orig = new std::ostringstream("");
		_descriptionLines[ section ] = oss_orig;
	}

	oss = new std::ostringstream("");
	addSpaces_( oss, parameterIndent_ );
	*oss << param;
	bool indentFirst = false;
	unsigned int charsUsed = oss->str().size();
	if (charsUsed > firstcolumnwidth) {
		*oss << std::endl;
		indentFirst = true;
	} else {
		addSpaces_( oss, firstcolumnwidth - charsUsed + 1 );
	}

	// Now we format the description, breaking it into lines
	std::vector< std::string > sublines;
	formatText_( sublines, description, firstcolumnwidth + 1, 0, maxWidth_ );
	if (!indentFirst) {
		sublines[ 0 ] = sublines[ 0 ].substr( firstcolumnwidth + 1 );
	}
	for (std::vector< std::string >::const_iterator it = sublines.cbegin();
		it != sublines.cend(); ++it) {
		*oss << *it << std::endl;
	}

	// Finally save the section names in order
	if (std::find( _sections.begin(), _sections.end(), section ) == _sections.end() ) {
		std::string tmpStr = section;
		_sections.push_back( tmpStr );
	}
	*oss_orig << oss->str();
	delete oss;
}

void NCPA::ParameterSet::addSpaces_( std::ostringstream *oss, unsigned int spaces ) {
	for (unsigned int i = 0; i < spaces; i++) {
		*oss << " ";
	}
}


/*
void NCPA::ParameterSet::addUsageLine( const std::string& line ) {
	std::string nline( line );
	_usage.push_back( nline );
}
*/

void NCPA::ParameterSet::printUsage( std::ostream& os ) const {

	// first the header
	std::vector< std::string >::const_iterator it;
	for (it = _headerLines.cbegin(); it != _headerLines.cend(); ++it ) {
		os << *it << std::endl;
	}
	os << std::endl;


	// now sections of parameters
	for (it = _sections.cbegin(); it != _sections.cend(); ++it) {
		os << *it << ":" << std::endl;
		std::ostringstream *oss = _descriptionLines.at( *it );
		os << oss->str()<< std::endl;

	}
	os << std::endl;

	// now the footer
	for (it = _footerLines.cbegin(); it != _footerLines.cend(); ++it ) {
		os << *it << std::endl;
	}
	os << std::endl;

/*
	for ( std::vector< std::string >::const_iterator it = _usage.begin();
		it != _usage.end(); ++it) {
		os << *it << std::endl;
	}
*/
}

void NCPA::ParameterSet::setDelimiters( std::string newdelim ) {
	_delims = newdelim;
}

void NCPA::ParameterSet::setComments( std::string newcomms ) {
	_comments = newcomms;
}

void NCPA::ParameterSet::addParameter( GenericParameter *newParam ) {
	_params.push_back( newParam );
}

void NCPA::ParameterSet::setStrict( bool tf ) {
	_strict = tf;
}

bool NCPA::ParameterSet::validate() {
	bool allgood = true;
	_failed_tests.clear();
	
	for ( std::vector< NCPA::ParameterTest * >::iterator it = _tests.begin(); 
			it != _tests.end(); ++it ) {
		bool passed = (*it)->validate( _params );
		if (! passed) {
			_failed_tests.push_back( *it );
		}
		//NCPA::GenericParameter *param = _params.findParameter( (*it)->optionName() );
		allgood = passed && allgood;
	}
	
	return allgood;
}

void NCPA::ParameterSet::printFailedTests( std::ostream& os ) const {
	for ( std::vector< NCPA::ParameterTest * >::const_iterator it = _failed_tests.begin(); 
			it != _failed_tests.end(); ++it ) {
		os << (*it)->failureMessage() << std::endl;
	}
}

unsigned int NCPA::ParameterSet::parseCommandLine( unsigned int argc, char **argv ) {
	//bool expectingArg = false;
	std::string lastarg;
	unsigned int nOptions = 0;
	
	for (unsigned int i = 1; i < argc; i++) {
		
		std::string currentarg = argv[ i ];
		if (isLongOption_( currentarg )) {
			i = processLongOption_( argc, argv, i );
			nOptions++;
		} else if (isShortOption_( currentarg )) {
			i = processShortOption_( argc, argv, i );
			nOptions++;
		} else {
			_unparsed.push_back( currentarg );
		}
	}
	return nOptions;
}

// returns true if the string is at least 3 characters and the first two are "--"
bool NCPA::ParameterSet::isLongOption_( std::string opt ) const {
	if (opt.size() < 3) {
		return false;
	}
	if (opt.compare(0,2,"--") == 0) {
		return true;
	} else {
		return false;
	}
}

// returns true if the string is at least 2 characters, the first character is "-", and the
// second is NOT "-"
bool NCPA::ParameterSet::isShortOption_( std::string opt ) const {
	if (opt.size() < 2) {
		return false;
	}
	if (opt[0] == '-') {
		if (opt[ 1 ] == '-') {
			return false;    // it's a long option
		}
		if (isdigit( opt[ 1 ] ) ) {
			return false;    // a negative number
		}
		return true;
	}
	return false;
}

unsigned int NCPA::ParameterSet::processShortOption_( int argc, char **argv,
	unsigned int i) {

	std::string fullarg = argv[ i ];

	// first, strip off the leading hyphen
	std::string stripped = fullarg.substr( fullarg.find_first_not_of( "-" ) );

	NCPA::GenericParameter *param = NULL;

	// for short options, we treat them each as a flag
	for (unsigned int i = 0; i < stripped.length(); i++) {
		std::string ch( 1, stripped[ i ] );
		param = _params.findParameter( ch );
		if (param != NULL) {
			param->setFound( true );
		} else {
			if (_strict) {
				throw std::invalid_argument( "Unexpected argument '-" + ch + "'." );
			} else {
				param = new FlagParameter( ch );
				param->setFound( true );
				_params.push_back( param );
			}
		}
	}
	return i;
}

void NCPA::ParameterSet::processDoubleOption_( std::string key, std::string value ) {

	// see of we're expecting this parameter
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param != NULL) {
		if (! param->needsArgument()) {
			throw std::invalid_argument( "Option " + key + " appears to have an argument '"
				+ value + "', but is not expecting one." );
		}
		param->parseArgument( value );
		param->setFound( true );
	} else {
		//not expecting it.  If in strict mode, throw an exception, otherwise
		// make a StringParameter
		if (_strict) {
			throw std::invalid_argument( "Unexpected argument: '" + key + "'" );
		} else {
			param = new StringParameter( key, value );
			param->setFound( true );
			_params.push_back( param );
		}
	}
}

void NCPA::ParameterSet::processSingleOption_( std::string key ) {

	// see if we're expecting this parameter
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param != NULL) {
		// expected it, let's see if it's expecting a value
		if (param->needsArgument()) {
			throw std::invalid_argument( "Option " + key + " expects an argument, but "
				+ "none appears to be provided." );
		}
		param->setFound( true );
	} else {
		// not expected.  strict mode?
		if (_strict) {
			throw std::invalid_argument( "Unexpected argument: '" + key + "'" );
		} else {
			param = new FlagParameter( key );
			param->setFound( true );
			_params.push_back( param );
		}
	}
}

// processes options in the form --flag, --option value, or --option=value, where the '='
// stands for any of the set delimiter characters
unsigned int NCPA::ParameterSet::processLongOption_( int argc, char **argv, 
	unsigned int i ) {
	
	std::string fullarg = argv[ i ];
	
	// first, strip off the leading hyphens
	std::string stripped = fullarg.substr( fullarg.find_first_not_of( "-" ) );
	
	// Two cases.  One argument with a delimiter character, or two arguments with spaces
	std::size_t found = stripped.find_first_of(_delims);
	std::string key, value;
	NCPA::GenericParameter *param = NULL;

	if ( found != std::string::npos ) {    // found one
		key = stripped.substr( 0, found );
		value = stripped.substr( found+1, std::string::npos );
		this->processDoubleOption_( key, value );
	} else {
		key = stripped;
		param = _params.findParameter( key );
		if (param != NULL) {
			
			// we're expecting it, now to see if it should have an argument with it
			if (param->needsArgument()) {
				
				// see if we can get the next value
				if ((int)i < (argc-1)) {
					// there are still arguments to parse out after this one
					value = argv[ i+1 ];

					// make sure the value isn't structured like another option or flag
					if (isLongOption_(value) || isShortOption_(value)) {
						std::ostringstream oss;
						oss << "Value '" << value << "' for option '" << key
							<< "' appears to be another argument.";
						throw std::invalid_argument( oss.str() );
					}

					// value looks legit, so use it
					this->processDoubleOption_( key, value );
					//param->parseArgument( value );
					//param->setFound( true );
					i++;

				} else {
					// expecting another argument but there aren't any more
					std::ostringstream oss;
					oss << "Option '" << key << "' expects a value, but none was provided.";
					throw std::invalid_argument( oss.str() );
				}
			} else {
				// no argument expected, so mark that it's here and move on
				this->processSingleOption_( key );
				//param->setFound( true );
			}
		} else {
			// We weren't expect this one.  If in strict mode, throw an exception
			if (_strict) {
				throw std::invalid_argument( "Unknown option '" + fullarg + "'." );
			} else {
				// don't know what to do with it, so let's see if it's got an argument
				if ((int)i < argc-1) {
					value = argv[ i+1 ];
					if (isLongOption_(value) || isShortOption_(value)) {
						// next option looks like an option, so treat this one as a flag
						this->processSingleOption_( key );
						//param = new FlagOption( key, true );
						//param->setFound( true );
						//_params.push_back( param );
					} else {
						// next option looks like a value, so use it as a string
						this->processDoubleOption_( key, value );
						//param = new StringOption( key, value );
						//param->setFound( true );
						//_params.push_back( param );
						i++;
					}
				} else {
					// no more arguments to parse, so it must be a flag
					this->processSingleOption_( key );
					//param = new FlagOption( key, true );
					//param->setFound( true );
					//_params.push_back( param );
				}		
			}
		}
	}

	return i;
}

unsigned int NCPA::ParameterSet::parseFile( std::string filename ) {

	std::ifstream ifs( filename, std::ifstream::in );
	if (! ifs.good() ) {
		return 0;
	}
	std::string line;
	//NCPA::GenericParameter *param = NULL;
	unsigned int linesParsed = 0;
	char linebuffer[ 1024 ];
	//memset(linebuffer,0,1024);

	do {
		std::getline( ifs, line );

		// see if it's blank
		line = NCPA::deblank( line );
		if (line.length() > 0 && line.find_first_of( _comments ) > 0) {

			// look for delimiters using strtok
			std::vector< std::string > tokens;

			char *pch;
			memset(linebuffer,0,1024);
			std::strcpy( linebuffer, line.c_str() );
			pch = strtok( linebuffer, _delims.c_str() );
			while (pch != NULL) {
				tokens.push_back( pch );
				pch = strtok( NULL, _delims.c_str() );
			}

			if (tokens.size() > 2) {
				throw std::invalid_argument( "Ambiguous line in file " + filename 
					+ ": '" + line + "'." );
			}

			if (tokens.size() == 1) {
				// no argument provided, just the flag
				this->processSingleOption_( NCPA::deblank(tokens[ 0 ]) );
			} else {
				this->processDoubleOption_( NCPA::deblank(tokens[ 0 ]),
				NCPA::deblank( tokens[ 1 ] ) );
			}

			linesParsed++;
		}
		line.clear();
	} while (ifs.good());

	ifs.close();

	return linesParsed;
}

NCPA::GenericParameter * NCPA::ParameterSet::getParameter( std::string key ) {
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param == NULL) {
		param = new NCPA::NullParameter( key );
		_params.push_back( param );
	}
	return param;
}


int NCPA::ParameterSet::getInteger( std::string key ) const {
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param == NULL) {
		throw std::invalid_argument( "No parameter '" + key 
			+ "' has been specified!" );
	}

	return param->getIntegerValue();
}

double NCPA::ParameterSet::getFloat( std::string key ) const {
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param == NULL) {
		throw std::invalid_argument( "No parameter '" + key 
			+ "' has been specified!" );
	}

	return param->getFloatValue();
}

std::string NCPA::ParameterSet::getString( std::string key ) const {
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param == NULL) {
		throw std::invalid_argument( "No parameter '" + key 
			+ "' has been specified!" );
	}

	return param->getStringValue();
}

bool NCPA::ParameterSet::getBool( std::string key ) const {
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param == NULL) {
		throw std::invalid_argument( "No parameter '" + key 
			+ "' has been specified!" );
	}

	return param->getBoolValue();
}


bool NCPA::ParameterSet::wasFound( std::string key ) const {
	NCPA::GenericParameter *param = _params.findParameter( key );
	if (param == NULL) {
		throw std::invalid_argument( "No parameter '" + key 
			+ "' has been specified!" );
	}

	return param->wasFound();
}


void NCPA::ParameterSet::removeParameter( std::string key ) {
	for ( std::vector< NCPA::GenericParameter * >::iterator it
		= _params.begin(); it != _params.end(); ++it ) {

		if ((*it)->getKey() == key) {
			_params.erase( it );
		}
	}
}


NCPA::ParameterTest * NCPA::ParameterSet::addTest( const std::string& option,
		NCPA::PARAMETER_TEST_TYPE test_type ) {
			
	NCPA::ParameterTest *crit;
	switch (test_type) {
		case PARAMETER_TEST_REQUIRED:
			crit = new NCPA::RequiredTest( option );
			break;
		case PARAMETER_TEST_REQUIRED_IF:
			crit = new NCPA::RequiredIfOtherIsPresentTest( option );
			break;
		case PARAMETER_TEST_RADIO_BUTTON:
			crit = new NCPA::RadioButtonTest( option );
			break;
		case PARAMETER_TEST_STRING_SET:
			crit = new NCPA::StringSetTest( option );
			break;
		case PARAMETER_TEST_INTEGER_POSITIVE:
			crit = new NCPA::IntegerGreaterThanTest( option );
			crit->addIntegerParameter( 0 );
			break;
		case PARAMETER_TEST_INTEGER_NEGATIVE:
			crit = new NCPA::IntegerLessThanTest( option );
			crit->addIntegerParameter( 0 );
			break;
		case PARAMETER_TEST_INTEGER_ZERO:
			crit = new NCPA::IntegerEqualToTest( option );
			crit->addIntegerParameter( 0 );
			break;
		case PARAMETER_TEST_INTEGER_NONZERO:
			crit = new NCPA::IntegerNotEqualToTest( option );
			crit->addIntegerParameter( 0 );
			break;
		case PARAMETER_TEST_FLOAT_POSITIVE:
			crit = new NCPA::FloatGreaterThanTest( option );
			crit->addFloatParameter( 0.0 );
			break;
		case PARAMETER_TEST_FLOAT_NEGATIVE:
			crit = new NCPA::FloatLessThanTest( option );
			crit->addFloatParameter( 0.0 );
			break;
		case PARAMETER_TEST_FLOAT_ZERO:
			crit = new NCPA::FloatEqualToTest( option );
			crit->addFloatParameter( 0.0 );
			break;
		case PARAMETER_TEST_FLOAT_NONZERO:
			crit = new NCPA::FloatNotEqualToTest( option );
			crit->addFloatParameter( 0.0 );
			break;
		case PARAMETER_TEST_INTEGER_GREATER_THAN:
			crit = new NCPA::IntegerGreaterThanTest( option );
			break;
		case PARAMETER_TEST_INTEGER_GREATER_THAN_OR_EQUAL:
			crit = new NCPA::IntegerGreaterThanOrEqualToTest( option );
			break;
		case PARAMETER_TEST_INTEGER_LESS_THAN:
			crit = new NCPA::IntegerLessThanTest( option );
			break;
		case PARAMETER_TEST_INTEGER_LESS_THAN_OR_EQUAL:
			crit = new NCPA::IntegerLessThanOrEqualToTest( option );
			break;
		case PARAMETER_TEST_INTEGER_EQUAL:
			crit = new NCPA::IntegerEqualToTest( option );
			break;
		case PARAMETER_TEST_INTEGER_NOT_EQUAL:
			crit = new NCPA::IntegerNotEqualToTest( option );
			break;
		case PARAMETER_TEST_FLOAT_GREATER_THAN:
			crit = new NCPA::FloatGreaterThanTest( option );
			break;
		case PARAMETER_TEST_FLOAT_GREATER_THAN_OR_EQUAL:
			crit = new NCPA::FloatGreaterThanOrEqualToTest( option );
			break;
		case PARAMETER_TEST_FLOAT_LESS_THAN:
			crit = new NCPA::FloatLessThanTest( option );
			break;
		case PARAMETER_TEST_FLOAT_LESS_THAN_OR_EQUAL:
			crit = new NCPA::FloatLessThanOrEqualToTest( option );
			break;
		case PARAMETER_TEST_FLOAT_EQUAL:
			crit = new NCPA::FloatEqualToTest( option );
			break;
		case PARAMETER_TEST_FLOAT_NOT_EQUAL:
			crit = new NCPA::FloatNotEqualToTest( option );
			break;
		case PARAMETER_TEST_STRING_MINIMUM_LENGTH:
			crit = new NCPA::StringMinimumLengthTest( option );
			break;
		case PARAMETER_TEST_STRING_MAXIMUM_LENGTH:
			crit = new NCPA::StringMaximumLengthTest( option );
			break;
		default:
			throw std::invalid_argument( "Undefined test requested" );
	}
	_tests.push_back( crit );
	return crit;
}

NCPA::ParameterTest *NCPA::ParameterSet::addTest( NCPA::ParameterTest *newTest ) {
	_tests.push_back( newTest );
	return newTest;
}

void NCPA::ParameterSet::printParameters( bool printTests, std::ostream& os ) const {
	for ( ParameterVector::const_iterator it = _params.begin();
		it != _params.end(); ++it) {
		os << (*it)->status() << std::endl;
		if (printTests) {
			std::vector< NCPA::ParameterTest * > subset = this->getTests( (*it)->getKey() );
			for ( std::vector< NCPA::ParameterTest * >::const_iterator it2 = subset.begin();
				it2 != subset.end(); ++it2) {
				os << "   TEST: " << (*it2)->description() << std::endl;
			}
		}
	}
}

std::vector< std::string > NCPA::ParameterSet::getUnparsedOptions() const {
	std::vector< std::string > u( _unparsed );
	return u;
}

std::vector< NCPA::ParameterTest * > NCPA::ParameterSet::getTests( std::string key ) const {
	std::vector< NCPA::ParameterTest * > subset;
	for (std::vector< NCPA::ParameterTest * >::const_iterator it = _tests.begin();
		it != _tests.end(); ++it) {
		if ( (*it)->optionName() == key ) {
			subset.push_back( *it );
		}
	}
	return subset;
}


/*
*********************************************************************************
Code for ParameterVector class
*********************************************************************************
*/
NCPA::GenericParameter * NCPA::ParameterVector::findParameter( const std::string& key ) const {
	NCPA::GenericParameter *param = NULL;
	for ( std::vector< NCPA::GenericParameter * >::const_iterator it = this->begin();
		it != this->end(); ++it ) {

		if ((*it)->getKey() == key) {
			return *it;
		}
	}
	return param;
}


/*
*********************************************************************************
Code for GenericParameter abstract base class
*********************************************************************************
*/
NCPA::GenericParameter::~GenericParameter() { }

void NCPA::GenericParameter::setFound( bool tf ) {
	_found = tf;
}

bool NCPA::GenericParameter::wasFound() const {
	return _found;
}

bool NCPA::GenericParameter::isValid() const {
	return _valid;
}

std::string NCPA::GenericParameter::getKey() const {
	return _key;
}

bool NCPA::GenericParameter::isNull() const {
	return false;
}

/*
*********************************************************************************
Code for NullParameter class
*********************************************************************************
*/
NCPA::NullParameter::NullParameter( std::string key ) {
	_key = key;
	_found = false;
	_valid = false;
}
bool NCPA::NullParameter::needsArgument() const {
	return false;
}
void NCPA::NullParameter::parseArgument( const std::string& arg ) {
	throw std::invalid_argument( "Parsing an argument for the null parameter is invalid" );
}
std::string NCPA::NullParameter::description() const {
	return "";
}
std::string NCPA::NullParameter::status() const {
	return "";
}
int NCPA::NullParameter::getIntegerValue() const {
	throw std::invalid_argument( "Parameter " + _key + " not set or defined" );
}
double NCPA::NullParameter::getFloatValue() const {
	throw std::invalid_argument( "Parameter " + _key + " not set or defined" );
}
std::string NCPA::NullParameter::getStringValue() const {
	throw std::invalid_argument( "Parameter " + _key + " not set or defined" );
}
bool NCPA::NullParameter::getBoolValue() const {
	throw std::invalid_argument( "Parameter " + _key + " not set or defined" );
}

bool NCPA::NullParameter::isNull() const {
	return true;
}

/*
*********************************************************************************
Code for IntegerParameter class
*********************************************************************************
*/
NCPA::IntegerParameter::IntegerParameter( std::string key ) {
	_key = key;
	_found = false;
	_valid = false;
	_value = 0;
}

NCPA::IntegerParameter::IntegerParameter( std::string key, int defaultValue ) {
	_key = key;
	_found = false;
	setValue( defaultValue );
}

void NCPA::IntegerParameter::setValue( int newval ) {
	_value = newval;
	_valid = true;
}

bool NCPA::IntegerParameter::needsArgument() const { 
	return true;
}

int NCPA::IntegerParameter::getIntegerValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _value;
}

double NCPA::IntegerParameter::getFloatValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return (double)_value;
}

std::string NCPA::IntegerParameter::getStringValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return std::to_string( _value );
}

bool NCPA::IntegerParameter::getBoolValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _value != 0;
}

void NCPA::IntegerParameter::parseArgument( const std::string& arg ) {
	_value = std::stoi( arg );
	_valid = true;
}

std::string NCPA::IntegerParameter::description() const {
	return _key + ": Integer";
}

std::string NCPA::IntegerParameter::status() const {
	std::ostringstream os;

	os << this->description();
		if ( this->isValid() ) {
			os << ", value = " << std::to_string( _value )
				<< ( this->wasFound() ? " [provided]" : " [default]" );
		} else {
			os << ", no value set";
		}
		return os.str();
}

/*
*********************************************************************************
Code for FloatParameter class
*********************************************************************************
*/
NCPA::FloatParameter::FloatParameter( std::string key ) {
	_key = key;
	_found = false;
	_valid = false;
	_value = 0.0;
}

NCPA::FloatParameter::FloatParameter( std::string key, double defaultValue ) {
	_key = key;
	_found = false;
	setValue( defaultValue );
}

void NCPA::FloatParameter::setValue( double newval ) {
	_value = newval;
	_valid = true;
}

bool NCPA::FloatParameter::needsArgument() const { 
	return true;
}

int NCPA::FloatParameter::getIntegerValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return (int)std::round(_value);
}

double NCPA::FloatParameter::getFloatValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _value;
}

std::string NCPA::FloatParameter::getStringValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return std::to_string( _value );
}

bool NCPA::FloatParameter::getBoolValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _value != 0.0;
}

void NCPA::FloatParameter::parseArgument( const std::string& arg ) {
	_value = std::stof( arg );
	_valid = true;
}

std::string NCPA::FloatParameter::description() const {
	return _key + ": Floating-Point Number";
}

std::string NCPA::FloatParameter::status() const {
	std::ostringstream os;

	os << this->description();
		if ( this->isValid() ) {
			os << ", value = " << std::to_string( _value )
				<< ( this->wasFound() ? " [provided]" : " [default]" );
		} else {
			os << ", no value set";
		}
		return os.str();
}

/*
*********************************************************************************
Code for StringParameter class
*********************************************************************************
*/
NCPA::StringParameter::StringParameter( std::string key ) {
	_key = key;
	_found = false;
	_valid = false;
	_value = "";
}

NCPA::StringParameter::StringParameter( std::string key, std::string defaultValue ) {
	_key = key;
	_found = false;
	setValue( defaultValue );
}

void NCPA::StringParameter::setValue( std::string newval ) {
	_value = newval;
	_valid = true;
}

bool NCPA::StringParameter::needsArgument() const { 
	return true;
}

int NCPA::StringParameter::getIntegerValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return std::stoi( _value );
}

double NCPA::StringParameter::getFloatValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return std::stof( _value );
}

std::string NCPA::StringParameter::getStringValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	std::string valcopy = _value;
	return valcopy;
}

bool NCPA::StringParameter::getBoolValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _value.size() > 0;
}

void NCPA::StringParameter::parseArgument( const std::string& arg ) {
	_value = arg;
	_valid = true;
}

std::string NCPA::StringParameter::description() const {
	return _key + ": String";
}

std::string NCPA::StringParameter::status() const {
	std::ostringstream os;

	os << this->description();
		if ( this->isValid() ) {
			os << ", value = " << _value
				<< ( this->wasFound() ? " [provided]" : " [default]" );
		} else {
			os << ", no value set";
		}
		return os.str();
}

/*
*********************************************************************************
Code for FlagParameter class
*********************************************************************************
*/
NCPA::FlagParameter::FlagParameter( std::string key ) {
	_key = key;
	_found = false;
	_valid = true;
}

NCPA::FlagParameter::FlagParameter( std::string key, bool defaultValue ) {
	_key = key;
	_found = defaultValue;
	_valid = true;
}

void NCPA::FlagParameter::setValue( bool newval ) {
	_value = newval;
}

bool NCPA::FlagParameter::needsArgument() const { 
	return false;
}

int NCPA::FlagParameter::getIntegerValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _found ? 1 : 0;
}

double NCPA::FlagParameter::getFloatValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _found ? 1.0 : 0.0;
}

std::string NCPA::FlagParameter::getStringValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _found ? "true" : "false";
}

bool NCPA::FlagParameter::getBoolValue() const {
	if ( ! _valid ) {
		throw std::invalid_argument( "Parameter " + _key + " not set" );
	}
	return _found;
}

void NCPA::FlagParameter::parseArgument( const std::string& arg ) {
	_found = true;
}

bool NCPA::FlagParameter::isValid() const {
	return _valid;
}

std::string NCPA::FlagParameter::description() const {
	return _key + ": Boolean Flag";
}

std::string NCPA::FlagParameter::status() const {
	std::ostringstream os;

	os << this->description();
		if ( this->wasFound() ) {
			os << ", flag set";
		} else {
			os << ", flag not set";
		}
	return os.str();
}

/*
*********************************************************************************
Code for ParameterTest abstract base class
*********************************************************************************
*/

NCPA::ParameterTest::ParameterTest() : _optName{ "" }, _ready{ false } {}
NCPA::ParameterTest::ParameterTest( std::string option_name, bool starts_ready ) 
	: _optName{ option_name }, _ready{ starts_ready } {}

// Destructor for ABC, must not be pure virtual
NCPA::ParameterTest::~ParameterTest() { }

std::string NCPA::ParameterTest::optionName() const {
	return _optName;
}

void NCPA::ParameterTest::addIntegerParameter( int param ) { }
void NCPA::ParameterTest::addFloatParameter( double param ) { }
void NCPA::ParameterTest::addStringParameter( std::string param ) { }
bool NCPA::ParameterTest::ready() const { return _ready; }



/**********************************************************************
NCPA::ParameterTest derived class methods here
Each gets a constructor, a validate() method, a description, a failure message,
and optional add<Type>Parameter() method(s)
**********************************************************************/

NCPA::RequiredTest::RequiredTest( const std::string& optionName ) 
	: ParameterTest( optionName, true ) {}

std::string NCPA::RequiredTest::description() const {
	return _optName + " is present.";
}

std::string NCPA::RequiredTest::failureMessage() const {
	return _optName + " is not present.";
}

bool NCPA::RequiredTest::validate( const ParameterVector& paramVec )  {
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );
	// Just check to see if it's been provided
	if (param != NULL && param->wasFound()) {
		return true;
	} 
	return false;
}

std::string NCPA::RequiredTest::valueString() const { return ""; }



NCPA::RequiredIfOtherIsPresentTest::RequiredIfOtherIsPresentTest( 
	const std::string& optionName ) 
	: ParameterTest( optionName, false ) {}

NCPA::RequiredIfOtherIsPresentTest::RequiredIfOtherIsPresentTest(
	const std::string& optionName, const std::string& prereq ) 
	: ParameterTest( optionName, false ) {
	std::string new_prereq{ prereq };
	_prereqs.push_back( prereq );
	_ready = true;
}

NCPA::RequiredIfOtherIsPresentTest::RequiredIfOtherIsPresentTest(
	const std::string& optionName, unsigned int nPrereqs, std::string *prereqs ) 
	: ParameterTest( optionName, false ) {
	std::string pr;
	for (unsigned int i = 0; i < nPrereqs; i++) {
		pr = prereqs[ i ];
		_prereqs.push_back( pr );
	}
	_ready = true;
}

NCPA::RequiredIfOtherIsPresentTest::RequiredIfOtherIsPresentTest(
	const std::string& optionName, const std::vector< std::string > prereq_vector ) 
	: ParameterTest( optionName, false ), _prereqs{ prereq_vector } {
	//_prereqs = prereq_vector;
	_ready = true;
}

std::string NCPA::RequiredIfOtherIsPresentTest::description() const {
	return _optName + " is present if one of " + this->valueString() 
		+ " is also present.";
}
std::string NCPA::RequiredIfOtherIsPresentTest::failureMessage() const {
	return "One of " + this->valueString() + " is set, but " + _optName
		+ " is not set.";
}
std::string NCPA::RequiredIfOtherIsPresentTest::valueString() const {
	std::ostringstream oss;
	oss << "{ ";
	for (std::vector<std::string>::const_iterator it = _prereqs.begin();
			it != _prereqs.end(); ++it) {
		if (it != _prereqs.begin()) {
			oss << ", ";
		}
		oss << *it;
	}
	oss << " }";
	return oss.str();
}
bool NCPA::RequiredIfOtherIsPresentTest::validate( const ParameterVector& paramVec )  {
	
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	// see if at least one prereq is present
	bool prereqs_met = false;
	NCPA::GenericParameter *prereq_param = NULL;
	for (std::vector<std::string>::const_iterator it = _prereqs.begin();
			it != _prereqs.end(); ++it) {
		prereq_param = paramVec.findParameter( *it );
		if ( prereq_param != NULL && prereq_param->wasFound() ) {
			prereqs_met = true;
		}
	}

	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	// if prereq(s) there, check for option presence, otherwise return true
	if (prereqs_met) {
		return (param != NULL && param->wasFound());
	} else {
		return true;
	}
}
void NCPA::RequiredIfOtherIsPresentTest::addStringParameter( const std::string param ) {
	std::string str = param;
	_prereqs.push_back( str );
}
bool NCPA::RequiredIfOtherIsPresentTest::ready() const {
	return !( _prereqs.empty() );
}



NCPA::RadioButtonTest::RadioButtonTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ) {
	_buttons.clear();
	_matched.clear();
}
NCPA::RadioButtonTest::RadioButtonTest( const std::string& optionName, 
	unsigned int nButtons, std::string *buttons ) 
	: ParameterTest( optionName, false ) {
	_buttons.clear();
	_matched.clear();
	for (unsigned int i = 0; i < nButtons; i++) {
		std::string but = buttons[ i ];
		_buttons.push_back( but );
	}
}
NCPA::RadioButtonTest::RadioButtonTest( const std::string& optionName,
	const std::vector< std::string > newButtons ) 
	: ParameterTest( optionName, false ), _buttons{ newButtons } {
	//_buttons.clear();
	//_buttons = newButtons;
	//_matched.clear();
}
std::string NCPA::RadioButtonTest::description() const {
	return _optName + ": One and only one of " + this->valueString() + " must be present.";
}
std::string NCPA::RadioButtonTest::failureMessage() const {
	std::ostringstream oss;
	oss << _optName << ": " << _matched.size() << " of " << this->valueString()
		<< " are present; must be one and only one.";
	return oss.str();
}
std::string NCPA::RadioButtonTest::valueString() const {
	std::ostringstream oss;
	oss << "{ ";
	for (std::vector<std::string>::const_iterator it = _buttons.begin();
			it != _buttons.end(); ++it) {
		if (it != _buttons.begin()) {
			oss << ", ";
		}
		oss << *it;
	}
	oss << " }";
	return oss.str();
}
bool NCPA::RadioButtonTest::validate( const ParameterVector& paramVec )  {
	_matched.clear();
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = NULL;
	for (std::vector<std::string>::const_iterator it = _buttons.begin();
			it != _buttons.end(); ++it) {
		param = paramVec.findParameter( *it );
		if ( param != NULL && param->wasFound() ) {
			_matched.push_back( *it );
		}
	}
	
	return (_matched.size() == 1);
}
void NCPA::RadioButtonTest::addStringParameter( const std::string newButton ) {
	std::string str = newButton;
	_buttons.push_back( str );
}
std::vector< std::string > NCPA::RadioButtonTest::lastMatched() const {
	std::vector< std::string > v( _matched );
	return v;
}
bool NCPA::RadioButtonTest::ready() const {
	return !( _buttons.empty() );
}





NCPA::IntegerGreaterThanTest::IntegerGreaterThanTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 }, _testedValue{ 0 } {}

NCPA::IntegerGreaterThanTest::IntegerGreaterThanTest( const std::string& optionName,
	int comparison )
 	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0 } {}

std::string NCPA::IntegerGreaterThanTest::description() const {
	return _optName + " is greater than " + 
		this->valueString() + " if present.";
}
std::string NCPA::IntegerGreaterThanTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be greater than " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerGreaterThanTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getIntegerValue();
	return (_testedValue > _value);
}
std::string NCPA::IntegerGreaterThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerGreaterThanTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerGreaterThanOrEqualToTest::IntegerGreaterThanOrEqualToTest( 
	const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 }, _testedValue{ 0 } {}

NCPA::IntegerGreaterThanOrEqualToTest::IntegerGreaterThanOrEqualToTest( 
	const std::string& optionName, int comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0 } {}

std::string NCPA::IntegerGreaterThanOrEqualToTest::description() const {
	return _optName + " is greater than or equal to " + 
		this->valueString() + " if present.";
}
std::string NCPA::IntegerGreaterThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) 
		+ ") must be greater than or equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerGreaterThanOrEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getIntegerValue();
	return (_testedValue >= _value);
}
std::string NCPA::IntegerGreaterThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerGreaterThanOrEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerLessThanTest::IntegerLessThanTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 }, _testedValue{ 0 } {}

NCPA::IntegerLessThanTest::IntegerLessThanTest( const std::string& optionName, int comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0 } {}

std::string NCPA::IntegerLessThanTest::description() const {
	return _optName + " is less than " + 
		this->valueString() + ".";
}
std::string NCPA::IntegerLessThanTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be less than " 
		+ this->valueString() + " if present.";
}
bool NCPA::IntegerLessThanTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getIntegerValue();
	return (_testedValue < _value);
}
std::string NCPA::IntegerLessThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerLessThanTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerLessThanOrEqualToTest::IntegerLessThanOrEqualToTest( 
	const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 }, _testedValue{ 0 } {}

NCPA::IntegerLessThanOrEqualToTest::IntegerLessThanOrEqualToTest( 
	const std::string& optionName, int comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0 } {}

std::string NCPA::IntegerLessThanOrEqualToTest::description() const {
	return _optName + " is less than or equal to " + 
		this->valueString() + " if present.";
}
std::string NCPA::IntegerLessThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) 
		+ ") must be less than or equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerLessThanOrEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getIntegerValue();
	return (_testedValue <= _value);
}
std::string NCPA::IntegerLessThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerLessThanOrEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}





NCPA::IntegerEqualToTest::IntegerEqualToTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 }, _testedValue{ 0 } {}

NCPA::IntegerEqualToTest::IntegerEqualToTest( 
	const std::string& optionName, int comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0 } {}

std::string NCPA::IntegerEqualToTest::description() const {
	return _optName + " is equal to " + this->valueString() + ".";
}
std::string NCPA::IntegerEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be equal to " 
		+ this->valueString() + " if present.";
}
bool NCPA::IntegerEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getIntegerValue();
	return (_testedValue == _value);
}
std::string NCPA::IntegerEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerNotEqualToTest::IntegerNotEqualToTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 }, _testedValue{ 0 } {}

NCPA::IntegerNotEqualToTest::IntegerNotEqualToTest( 
	const std::string& optionName, int comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0 } {}

std::string NCPA::IntegerNotEqualToTest::description() const {
	return _optName + " is not equal to " + this->valueString() + " if present.";
}
std::string NCPA::IntegerNotEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must not be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerNotEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getIntegerValue();
	return (_testedValue != _value);
}
std::string NCPA::IntegerNotEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerNotEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}






NCPA::FloatGreaterThanTest::FloatGreaterThanTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0.0 }, _testedValue{ 0.0 } {}

NCPA::FloatGreaterThanTest::FloatGreaterThanTest( 
	const std::string& optionName, double comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0.0 } {}

std::string NCPA::FloatGreaterThanTest::description() const {
	return _optName + " is greater than " 
		 + this->valueString() + " if present.";
}
std::string NCPA::FloatGreaterThanTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be greater than " 
		+ this->valueString() + ".";
}
bool NCPA::FloatGreaterThanTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getFloatValue();
	return (_testedValue > _value);
}
std::string NCPA::FloatGreaterThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatGreaterThanTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}



NCPA::FloatGreaterThanOrEqualToTest::FloatGreaterThanOrEqualToTest( 
	const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0.0 }, _testedValue{ 0.0 } {}

NCPA::FloatGreaterThanOrEqualToTest::FloatGreaterThanOrEqualToTest( 
	const std::string& optionName, double comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0.0 } {}

std::string NCPA::FloatGreaterThanOrEqualToTest::description() const {
	return _optName + " is greater than or equal to " 
		 + this->valueString() + " if present.";
}
std::string NCPA::FloatGreaterThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) 
		+ ") must be greater than or equal to " 
		+ this->valueString() + ".";
}
bool NCPA::FloatGreaterThanOrEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getFloatValue();
	return (_testedValue >= _value);
}
std::string NCPA::FloatGreaterThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatGreaterThanOrEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}





NCPA::FloatLessThanTest::FloatLessThanTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0.0 }, _testedValue{ 0.0 } {}

NCPA::FloatLessThanTest::FloatLessThanTest( 
	const std::string& optionName, double comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0.0 } {}

std::string NCPA::FloatLessThanTest::description() const {
	return _optName + " is less than " 
		 + this->valueString() + " if present.";
}
std::string NCPA::FloatLessThanTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be less than " 
		+ this->valueString() + ".";
}
bool NCPA::FloatLessThanTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getFloatValue();
	return (_testedValue < _value);
}
std::string NCPA::FloatLessThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatLessThanTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}



NCPA::FloatLessThanOrEqualToTest::FloatLessThanOrEqualToTest( 
	const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0.0 }, _testedValue{ 0.0 } {}

NCPA::FloatLessThanOrEqualToTest::FloatLessThanOrEqualToTest( 
	const std::string& optionName, double comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0.0 } {}

std::string NCPA::FloatLessThanOrEqualToTest::description() const {
	return _optName + " is less than " 
		 + this->valueString() + " if present.";
}
std::string NCPA::FloatLessThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be less than " 
		+ this->valueString() + ".";
}
bool NCPA::FloatLessThanOrEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getFloatValue();
	return (_testedValue <= _value);
}
std::string NCPA::FloatLessThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatLessThanOrEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}




NCPA::FloatEqualToTest::FloatEqualToTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0.0 }, _testedValue{ 0.0 } {}

NCPA::FloatEqualToTest::FloatEqualToTest( 
	const std::string& optionName, double comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0.0 } {}

std::string NCPA::FloatEqualToTest::description() const {
	return _optName + " is equal to " + this->valueString() + " if present.";
}
std::string NCPA::FloatEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::FloatEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getFloatValue();
	return (_testedValue == _value);
}
std::string NCPA::FloatEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}



NCPA::FloatNotEqualToTest::FloatNotEqualToTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0.0 }, _testedValue{ 0.0 } {}

NCPA::FloatNotEqualToTest::FloatNotEqualToTest( 
	const std::string& optionName, double comparison ) 
	: ParameterTest( optionName, true ), _value{ comparison }, _testedValue{ 0.0 } {}

std::string NCPA::FloatNotEqualToTest::description() const {
	return _optName + " is not equal to " + this->valueString() + " if present.";
}
std::string NCPA::FloatNotEqualToTest::failureMessage() const {
	return _optName + " (" + std::to_string(_testedValue) + ") must not be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::FloatNotEqualToTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getFloatValue();
	return (_testedValue != _value);
}
std::string NCPA::FloatNotEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatNotEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}






NCPA::StringMinimumLengthTest::StringMinimumLengthTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 } {}

NCPA::StringMinimumLengthTest::StringMinimumLengthTest( const std::string& optionName,
	size_t minlength ) 
	: ParameterTest( optionName, true ), _value{ minlength } {}

std::string NCPA::StringMinimumLengthTest::description() const {
	return _optName + " is at least " + this->valueString() + " characters if present.";
}
std::string NCPA::StringMinimumLengthTest::failureMessage() const {
	return _optName + " (\"" + _testedValue + "\") must be at least " 
		+ this->valueString() + " characters long.";
}
bool NCPA::StringMinimumLengthTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getStringValue();
	return (_testedValue.size() >= _value);
}
std::string NCPA::StringMinimumLengthTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::StringMinimumLengthTest::addIntegerParameter( int param ) {
	if (param < 0) {
		throw std::range_error( "String length must not be negative" );
	}
	_value = (size_t)param;
	_ready = true;
}




NCPA::StringMaximumLengthTest::StringMaximumLengthTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ), _value{ 0 } {}

NCPA::StringMaximumLengthTest::StringMaximumLengthTest( const std::string& optionName,
	size_t maxlength ) 
	: ParameterTest( optionName, true ), _value{ maxlength } {}

std::string NCPA::StringMaximumLengthTest::description() const {
	return _optName + " is at most " + this->valueString() + " characters if present.";
}
std::string NCPA::StringMaximumLengthTest::failureMessage() const {
	return _optName + " (\"" + _testedValue + "\") must be at most " 
		+ this->valueString() + " characters long.";
}
bool NCPA::StringMaximumLengthTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getStringValue();
	return (_testedValue.size() <= _value);
}
std::string NCPA::StringMaximumLengthTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::StringMaximumLengthTest::addIntegerParameter( int param ) {
	if (param < 0) {
		throw std::range_error( "String length must not be negative" );
	}
	_value = (size_t)param;
	_ready = true;
}




NCPA::StringSetTest::StringSetTest( const std::string& optionName ) 
	: ParameterTest( optionName, false ) {}

NCPA::StringSetTest::StringSetTest( const std::string& optionName, unsigned int nSet,
	std::string *choices ) 
	: ParameterTest( optionName, false ) {
	//_choices.clear();
	for (unsigned int i = 0; i < nSet; i++) {
		std::string ch = choices[ i ];
		_choices.push_back( ch );
	}
}
NCPA::StringSetTest::StringSetTest( const std::string& optionName, 
	std::vector< std::string > choices ) 
	: ParameterTest( optionName, false ), _choices{ choices } {}

std::string NCPA::StringSetTest::description() const {
	return _optName + " must be in " + this->valueString() + " if present.";
}
std::string NCPA::StringSetTest::failureMessage() const {
	return _optName + ": " + '"' + _testedValue + '"' + " is not in " + this->valueString() + ".";
}
std::string NCPA::StringSetTest::valueString() const {
	std::ostringstream oss;
	oss << "{ ";
	for (std::vector<std::string>::const_iterator it = _choices.begin();
			it != _choices.end(); ++it) {
		if (it != _choices.begin()) {
			oss << ", ";
		}
		oss << '"' << *it << '"';
	}
	oss << " }";
	return oss.str();
}
bool NCPA::StringSetTest::validate( const ParameterVector& paramVec )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	NCPA::GenericParameter *param = paramVec.findParameter( _optName );

	if (param == NULL || ( ! param->isValid() ) ) {
		return true;
	}

	_testedValue = param->getStringValue();
	std::vector< std::string >::const_iterator it = std::find( 
		_choices.begin(), _choices.end(), _testedValue );
	return ( it != _choices.end() );
}
void NCPA::StringSetTest::addStringParameter( std::string newChoice ) {
	_choices.push_back( newChoice );
}
bool NCPA::StringSetTest::ready() const {
	return ! ( _choices.empty() );
}
