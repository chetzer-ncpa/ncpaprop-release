#ifndef NCPA_GEOGRAPHIC_H_INCLUDED
#define NCPA_GEOGRAPHIC_H_INCLUDED

#include "util.h"

/**
 * @version 1.2.1
 * @date 2021-08-03
 */

namespace NCPA {

	class Location {
		protected:
			double lat_, lon_, elev_;

		public:
			Location();
			Location( double lat, double lon, double elev = 0 );
			bool operator==( const Location &other ) const;
			bool operator!=( const Location &other ) const;
			bool operator<( const Location &other ) const;
			bool operator>( const Location &other ) const;
			double lat() const;
			double lon() const;
			double elev() const;
			void setLat( double newlat );
			void setLon( double newlon );
			void setElev( double newelev );
	};

	double azimuth( double lat1, double lon1, double lat2, double lon2 );
	double backazimuth( double lat1, double lon1, double lat2, double lon2 );
	void great_circle( double startlat, double startlon, double azimuth,
			   double range, int length, double *pathlat,
			   double *pathlon );
	double range( double lat1, double lon1, double lat2, double lon2 );
	void intersection( double lat1_deg, double lon1_deg, double bearing1_deg,
			   double lat2_deg, double lon2_deg, double bearing2_deg,
			   double &lat3_deg, double &lon3_deg );
	//double sphazimuth( double lat1, double lon1, double lat2, double lon2 );
	//double sphrange( double lat1, double lon1, double lat2, double lon2 );
	double earthradius( double lat );
	Location xy2latlon( double x, double y, double lat0, double lon0 );
	void xy2latlon( double x, double y, double lat0, double lon0, double &newlat, double &newlon );
	//Location xyz2latlon( double x, double y, double z, double r );
	void latlon2xy( size_t npts, double *lat, double *lon,
		double *x, double *y, double reflat = -999.9 );
	void latlon2xy( const std::vector<double> &lat, const std::vector<double> &lon,
		std::vector<double> &x, std::vector<double> &y, double reflat = -999.9 );

	//angular distance to linear
	double deg2km( double angular );
	double km2deg( double linear );
	
	std::string dms_str( double coord );
	double normalizeLon( double lon );
}


#endif
