#ifndef NCPAPROP_MESSAGING_H_INCLUDED
#define NCPAPROP_MESSAGING_H_INCLUDED

#include <iostream>

namespace NCPA {
	enum class VerbosityLevel {
		QUIET=0, ERROR, INFO, DEBUG
	};


	class Messenger {

		public:
			Messenger();
			Messenger(NCPA::VerbosityLevel v);
			Messenger(NCPA::VerbosityLevel v,
					std::ostream &errorstr );
			Messenger(NCPA::VerbosityLevel v,
					std::ostream &errorstr,
					std::ostream &infostr );
			Messenger(NCPA::VerbosityLevel v,
					std::ostream &errorstr,
					std::ostream &infostr,
					std::ostream &debugstr );

			void debug( const std::string &msg );
			void error( const std::string &msg );
			void info(  const std::string &msg );

			void set_verbosity( NCPA::VerbosityLevel v );
			NCPA::VerbosityLevel get_verbosity() const;

		protected:
			NCPA::VerbosityLevel _verbosity;
			std::ostream &_error_stream;
			std::ostream &_info_stream;
			std::ostream &_debug_stream;

			void _echo_message( const std::string &msg, std::ostream &ostr );
	};
}


#endif
