#include "messaging.h"
#include <iostream>



NCPA::Messager::Messager(
			NCPA::VerbosityLevel v,
			std::ostream &errorstr,
			std::ostream &infostr,
			std::ostream &debugstr ) :
		_verbosity{v},
		_error_stream{errorstr},
		_info_stream{infostr},
		_debug_stream{debugstr}
	{}


NCPA::Messager::Messager() : Messager(
		NCPA::VerbosityLevel::INFO,
		std::cerr,
		std::cout,
		std::cout )
	{}

NCPA::Messager::Messager(NCPA::VerbosityLevel v) : Messager(
		v,
		std::cerr,
		std::cout,
		std::cout )
	{}

NCPA::Messager::Messager(NCPA::VerbosityLevel v, std::ostream &estr) : Messager(
		v,
		estr,
		std::cout,
		std::cout )
	{}

NCPA::Messager::Messager(
		NCPA::VerbosityLevel v,
		std::ostream &estr,
		std::ostream &istr ) : Messager(
		v,
		estr,
		istr,
		std::cout )
	{}

void NCPA::Messager::_echo_message( const std::string &msg, std::ostream &ostr ) {
	ostr << msg << std::endl;
}

void NCPA::Messager::debug( const std::string &msg ) {
	if (_verbosity >= NCPA::VerbosityLevel::DEBUG) {
		this->_echo_message(msg, _debug_stream);
	}
}

void NCPA::Messager::info( const std::string &msg ) {
	if (_verbosity >= NCPA::VerbosityLevel::INFO) {
		this->_echo_message(msg, _info_stream);
	}
}

void NCPA::Messager::error( const std::string &msg ) {
	if (_verbosity >= NCPA::VerbosityLevel::ERROR) {
		this->_echo_message(msg, _error_stream);
	}
}
