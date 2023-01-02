#include "messaging.h"
#include <iostream>



NCPA::Messenger::Messenger(
			NCPA::VerbosityLevel v,
			std::ostream &errorstr,
			std::ostream &infostr,
			std::ostream &debugstr ) :
		_verbosity{v},
		_error_stream{errorstr},
		_info_stream{infostr},
		_debug_stream{debugstr}
	{}


NCPA::Messenger::Messenger() : Messenger(
		NCPA::VerbosityLevel::INFO,
		std::cerr,
		std::cout,
		std::cout )
	{}

NCPA::Messenger::Messenger(NCPA::VerbosityLevel v) : Messenger(
		v,
		std::cerr,
		std::cout,
		std::cout )
	{}

NCPA::Messenger::Messenger(NCPA::VerbosityLevel v, std::ostream &estr) : Messenger(
		v,
		estr,
		std::cout,
		std::cout )
	{}

NCPA::Messenger::Messenger(
		NCPA::VerbosityLevel v,
		std::ostream &estr,
		std::ostream &istr ) : Messenger(
		v,
		estr,
		istr,
		std::cout )
	{}

void NCPA::Messenger::_echo_message( const std::string &msg, std::ostream &ostr ) {
	ostr << msg << std::endl;
}

void NCPA::Messenger::debug( const std::string &msg ) {
	if (_verbosity >= NCPA::VerbosityLevel::DEBUG) {
		this->_echo_message(msg, _debug_stream);
	}
}

void NCPA::Messenger::info( const std::string &msg ) {
	if (_verbosity >= NCPA::VerbosityLevel::INFO) {
		this->_echo_message(msg, _info_stream);
	}
}

void NCPA::Messenger::error( const std::string &msg ) {
	if (_verbosity >= NCPA::VerbosityLevel::ERROR) {
		this->_echo_message(msg, _error_stream);
	}
}
