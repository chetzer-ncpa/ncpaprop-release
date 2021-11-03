#ifndef NCPAPROP_PROFILESERIESATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_PROFILESERIESATMOSPHERE2D_H_INCLUDED

#include "Atmosphere2D.h"


namespace NCPA {
	class ProfileSeriesAtmosphere2D : public Atmosphere2D {

	public:
		ProfileSeriesAtmosphere2D();
		ProfileSeriesAtmosphere2D( const std::string &filename );
		ProfileSeriesAtmosphere2D( const std::string &filename,
			const std::string &headerfilename );
		ProfileSeriesAtmosphere2D( const std::string &filename,
			size_t skiplines );
		ProfileSeriesAtmosphere2D( const std::string &filename,
			const std::string &headerfilename,
			size_t skiplines );
		~ProfileSeriesAtmosphere2D();

	protected:
		void read_summary_file( const std::string &filename,
			std::vector< std::string > &atmlines );
		void process_summary_file_lines(
			const std::vector< std::string > &atmlines,
			const std::string &headerfilename, size_t skiplines );
	};
}




#endif