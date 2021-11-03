#include "ProfileSeriesAtmosphere2D.h"
#include "Atmosphere1D.h"
#include "units.h"
#include <fstream>


NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D() : Atmosphere2D() { }

NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D(
	const std::string& filename ) : Atmosphere2D() {

	std::vector< std::string > atmlines;
	read_summary_file( filename, atmlines );
	process_summary_file_lines( atmlines, "", 0 );
}

NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D( const std::string& filename, 
	const std::string &headerfilename ) : Atmosphere2D() {

	std::vector< std::string > atmlines;
	read_summary_file( filename, atmlines );
	process_summary_file_lines( atmlines, headerfilename, 0 );
}

NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D(
	const std::string& filename,
	size_t skiplines ) : Atmosphere2D() {

	std::vector< std::string > atmlines;
	read_summary_file( filename, atmlines );
	process_summary_file_lines( atmlines, "", skiplines );
}

NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D(
	const std::string& filename,
	const std::string &headerfilename,
	size_t skiplines ) : Atmosphere2D() {

	std::vector< std::string > atmlines;
	read_summary_file( filename, atmlines );
	process_summary_file_lines( atmlines, headerfilename, skiplines );
}



void NCPA::ProfileSeriesAtmosphere2D::process_summary_file_lines(
		const std::vector< std::string > &atmlines,
		const std::string &headerfilename, size_t skiplines ) {

	double range;
	std::ostringstream oss;
	std::vector< std::string > fields;
	NCPA::Atmosphere1D *tempatm;
	std::string atmfile;

	this->set_insert_range_units( NCPA::Units::fromString( "km" ) );
	for (std::vector< std::string >::const_iterator it = atmlines.begin(); it != atmlines.end(); ++it) {
		fields = NCPA::split( NCPA::deblank( *it ), " ," );
		if ( fields.size() != 2 )  {
			oss << "ProfileSeriesAtmosphere2D - Error parsing input line:" << std::endl << *it << std::endl
				<< "Must be formatted as:" << std::endl
				<< "range  filename" << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		try {
			range = std::stof( fields[ 0 ] );
		} catch ( std::invalid_argument &e ) {
			oss << "ProfileSeriesAtmosphere2D - Error parsing input line:" << std::endl << *it << std::endl
				<< "First field not parseable as a floating-point number" << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		atmfile = fields[ 1 ];
		if (headerfilename.size() > 0) {
			tempatm = new Atmosphere1D( atmfile, headerfilename, skiplines );
		} else {
			tempatm = new Atmosphere1D( atmfile, skiplines );
		}
		insert_profile( tempatm, range );
	}
	sort_profiles();
}

void NCPA::ProfileSeriesAtmosphere2D::read_summary_file(
		const std::string &filename,
		std::vector< std::string > &atmlines ) {

	std::string line;
	std::ifstream in( filename );
	std::getline( in, line );
	while ( in.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = NCPA::deblank( line );
		if (line[ 0 ] != '#') {
			atmlines.push_back( line );
		}

		getline( in, line );
	}
	in.close();

}





	// std::ifstream infile( filename );
	// double range;
	// std::string atmfile;
	// infile >> range >> atmfile;
	// Atmosphere1D *tempatm;
	// set_insert_range_units( NCPA::Units::fromString( "km" ) );
	// while (infile.good()) {
	// 	tempatm = new Atmosphere1D( atmfile );
	// 	insert_profile( tempatm, range );
	// 	infile >> range >> atmfile;
	// }
	// infile.close();
	// sort_profiles();


NCPA::ProfileSeriesAtmosphere2D::~ProfileSeriesAtmosphere2D() { }

//NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D( const ProfileSeriesAtmosphere2D &atm ) : Atmosphere2D( atm ) { }