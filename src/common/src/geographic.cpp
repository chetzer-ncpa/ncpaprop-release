#include <cstdio>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include "geographic.h"
#include "util.h"

#ifndef PI
#define PI 3.14159
#endif

// Geographic functions and such.

/**
 * Some geographic calculations copyright 2002-2011 Chris Veness, 
 * offered under a Creative Commons attribution license.  Original
 * work found at http://www.movable-type.co.uk/scripts/latlong.html
 */


NCPA::Location::Location( double lat, double lon, double elev ) {
	lat_ = lat;
	lon_ = lon;
	elev_ = elev;
}

NCPA::Location::Location() {
	lat_ = -999;
	lon_ = -999;
	elev_ = 0;
}

void NCPA::Location::setLat( double newlat ) {
	if (std::fabs( newlat ) > 90.0) {
		throw std::invalid_argument( "Latitude must be in the range [-90,90]!");
	}
	lat_ = newlat;
}

void NCPA::Location::setLon( double newlon ) {
	if (std::fabs( newlon ) > 180.0) {
		throw std::invalid_argument( "Longitude must be in the range [-180,180]!");
	}
	lon_ = newlon;
}

void NCPA::Location::setElev( double newelev ) {
	elev_ = newelev;
}

double NCPA::Location::lat() const { return lat_; }
double NCPA::Location::lon() const { return lon_; }
double NCPA::Location::elev() const { return elev_; }

bool NCPA::Location::operator==( const NCPA::Location &other ) const {
	if (lat_ == other.lat() && lon_ == other.lon() 
			&& elev_ == other.elev())
		return true;
	return false;
}

bool NCPA::Location::operator!=( const NCPA::Location &other ) const {
	return !((*this) == other);
}

bool NCPA::Location::operator<( const NCPA::Location &other ) const {
	if (lat_ < other.lat())
		return true;
	else if (lat_ > other.lat())
		return false;
	else if (lon_ < other.lon())
		return true;
	else
		return false;
}

bool NCPA::Location::operator>( const NCPA::Location &other ) const {
	if (lat_ > other.lat())
		return true;
	else if (lat_ < other.lat())
		return false;
	else if (lon_ > other.lon())
		return true;
	else
		return false;
}


double NCPA::azimuth( double lat1_deg, double lon1_deg, double lat2_deg, double lon2_deg ) {

	double dlon = NCPA::deg2rad( lon2_deg - lon1_deg);
	double lat1 = NCPA::deg2rad( lat1_deg );
	double lat2 = NCPA::deg2rad( lat2_deg );
	double y = std::sin( dlon ) * std::cos( lat2 );
	double x = std::cos( lat1 ) * std::sin( lat2 ) - std::sin( lat1 ) * std::cos( lat2 ) * std::cos( dlon );
	double bearing = NCPA::rad2deg( std::atan2( y, x ) );
	return NCPA::normalizeAzimuth( bearing );
	
}

double NCPA::backazimuth( double lat1_deg, double lon1_deg, double lat2_deg, double lon2_deg ) {

	return NCPA::azimuth( lat2_deg, lon2_deg, lat1_deg, lon1_deg );

}

std::string NCPA::dms_str( double coord ) {
	
	int sign = coord >= 0 ? 1 : -1;
	coord = std::fabs( coord );
	int degs = (int)(std::floor( coord ));
	coord -= degs;
	coord *= 60.0;
	int mins = (int)(std::floor( coord ));
	coord -= mins;
	coord *= 60.0;
	degs *= sign;
	
	char buffer[ 256 ];
	std::sprintf( buffer, "%dd%02d'%06.3f\"", degs, mins, coord );
	std::string s( buffer );
	return s;
}

double NCPA::normalizeLon(double lon) {
	while (lon > 180.0) {
		lon -= 360.0;
	}
	while (lon <= -180.0) {
		lon += 360.0;
	}
	return lon;
}


double NCPA::deg2km( double angular ) {
	return 2 * PI * 6371 / 360.0 * angular;
}

double NCPA::km2deg( double linear ) {
	return linear * 360.0 / (2.0 * PI * 6371 );
}

void NCPA::great_circle( double startlat, double startlon, double azimuth, 
		double range_km, int length, double *pathlat, double *pathlon ) {
	
	double lat1 = NCPA::deg2rad( startlat );
	double lon1 = NCPA::deg2rad( startlon );
	//double d = range_km / 6371.0;
	double brng = NCPA::deg2rad( azimuth );
	pathlat[0] = startlat;
	pathlon[0] = startlon;
	for (int i = 1; i < length; i++) {
		double d = (range_km / (length-1)) * i / 6371.0;
		pathlat[ i ] = NCPA::rad2deg( 
				std::asin( std::sin(lat1) * std::cos(d) + std::cos(lat1) * std::sin(d) * std::cos(brng) ) 
			);
		pathlon[ i ] = NCPA::normalizeLon( NCPA::rad2deg( 
				lon1 + std::atan2( std::sin(brng) * std::sin(d) * std::cos(lat1),
						   std::cos(d) - std::sin(lat1) * std::sin(NCPA::deg2rad(pathlat[i])) )
				     ) );
	}
	
}

void NCPA::intersection(double lat1_deg, double lon1_deg, double bearing1_deg, double lat2_deg, double lon2_deg, 
			double bearing2_deg, double& lat3_deg, double& lon3_deg) {

	double lat1 = NCPA::deg2rad( lat1_deg );
	double lat2 = NCPA::deg2rad( lat2_deg );
	double lon1 = NCPA::deg2rad( lon1_deg );
	double lon2 = NCPA::deg2rad( lon2_deg );
	double theta1 = NCPA::deg2rad( bearing1_deg );
	double theta2 = NCPA::deg2rad( bearing2_deg );
	double dlat = lat2 - lat1;
	double dlon = lon2 - lon1;
	
	double d12 = 2.0 * std::asin( 
		std::sqrt( std::sin(dlat/2)*std::sin(dlat/2) + std::cos(lat1)*std::cos(lat2)*std::sin(dlon/2)*std::sin(dlon/2) ) 
			);
	double phi1 = std::acos( std::sin(lat2) - std::sin(lat1)*std::cos(d12) / std::sin(d12)*std::cos(lat1) );
	double phi2 = std::acos( std::sin(lat1) - std::sin(lat2)*std::cos(d12) / std::sin(d12)*std::cos(lat2) );
	
	double theta12, theta21;
	if (std::sin(lon2 - lon1) > 0) {
		theta12 = phi1;
		theta21 = 2*PI - phi2;
	} else {
		theta12 = 2*PI - phi1;
		theta21 = phi2;
	}
	
	double alpha1 = std::fmod( theta1 - theta12 + PI, 2*PI ) - PI;
	double alpha2 = std::fmod( theta21 - theta2 + PI, 2*PI ) - PI;
	double alpha3 = std::acos( -std::cos( alpha1 )*std::cos( alpha2 ) + std::sin( alpha1 )*std::sin( alpha2 )*std::cos( d12 ) );
	double d13 = std::atan2( std::sin(d12)*std::sin(alpha1)*std::sin(alpha2), std::cos(alpha2) + std::cos(alpha1)*std::cos(alpha3) );
	double lat3 = std::asin( std::sin(lat1)*std::cos(d13) + std::cos(lat1)*std::sin(d13)*std::cos(theta1) );
	lat3_deg = NCPA::rad2deg( lat3 );
	double dlon13 = std::atan2( std::sin(theta1)*std::sin(d13)*std::cos(lat1), std::cos(d13) - std::sin(lat1)*std::sin(lat3) );
	lon3_deg = NCPA::rad2deg( std::fmod( lon1 + dlon13 + PI, 2*PI ) - PI );
}



/*
 * Computes the range along the surface of the earth using the Haversine formula
 */
double NCPA::range( double lat1_deg, double lon1_deg, double lat2_deg, double lon2_deg ) {

	double dLat = NCPA::deg2rad( lat2_deg - lat1_deg );
	double dLon = NCPA::deg2rad( lon2_deg - lon1_deg );
	double lat1 = NCPA::deg2rad( lat1_deg );
	double lat2 = NCPA::deg2rad( lat2_deg );
	
	double a = std::sin( dLat / 2 ) * std::sin( dLat / 2 ) +
			std::sin( dLon / 2 ) * std::sin( dLon / 2 ) * std::cos( lat1 ) * std::cos( lat2 );
	double c = 2.0 * std::atan2( std::sqrt( a ), std::sqrt( 1-a ) );
	return c * NCPA::earthradius( 0.5*lat1_deg + 0.5*lat2_deg );
	
}

double NCPA::earthradius( double lat ) {
/*
	double a = 6378.137;
	double b = 6356.7523;
	lat = deg2rad( lat );
	double numerator = a*a*a*a * std::cos( lat ) * std::cos( lat )
			 + b*b*b*b * std::sin( lat ) * std::sin( lat );
	double denominator = a*a * std::cos( lat )
			   + b*b * std::sin( lat );
	return std::sqrt( numerator / denominator );
*/
	return 6371.009;
}

void NCPA::xy2latlon( double x, double y, double lat0, double lon0, double &newlat, double &newlon ) {
	
	// Simple: go north by y km, then east by x km
	double templat[2], templon[2];
	double northaz = 0;
	double eastaz = 90;
	if (y < 0) {
		northaz = 180;
		y = -y;
	}
	if (x < 0) {
		eastaz = 270;
		x = -x;
	}
	
	NCPA::great_circle( lat0, lon0, northaz, y, 2, templat, templon );
	double templat2 = templat[1];
	double templon2 = templon[1];
	NCPA::great_circle( templat2, templon2, eastaz, x, 2, templat, templon );
	newlat = templat[1];
	newlon = templon[1];
	//Location loc( templat[1], templon[1] );
	//return loc;
}
	

NCPA::Location NCPA::xy2latlon( double x, double y, double lat0, double lon0 ) {
	
	// Simple: go north by y km, then east by x km
	double templat[2], templon[2];
	double northaz = 0;
	double eastaz = 90;
	if (y < 0) {
		northaz = 180;
		y = -y;
	}
	if (x < 0) {
		eastaz = 270;
		x = -x;
	}
	
	NCPA::great_circle( lat0, lon0, northaz, y, 2, templat, templon );
	double templat2 = templat[1];
	double templon2 = templon[1];
	NCPA::great_circle( templat2, templon2, eastaz, x, 2, templat, templon );
	Location loc( templat[1], templon[1] );
	return loc;
	
}

// convert a series of lat/lon points to absolute x/y points using an
// equirectangular projection, with reference latitide chosen to be the
// center of the latitude range
void NCPA::latlon2xy( size_t npts, double *lat, double *lon,
		double *x, double *y, double reflat ) {

	if (reflat < -90.0) {
		reflat = 0.5 * ( NCPA::max( lat, npts ) - NCPA::min( lat, npts ) );
	}
	double Re = NCPA::earthradius( reflat );
	for (size_t i = 0; i < npts; i++) {
		x[ i ] = Re * NCPA::deg2rad( lon[ i ] ) * std::cos( NCPA::deg2rad( reflat ) );
		y[ i ] = Re * NCPA::deg2rad( lat[ i ] );
	}
}

void NCPA::latlon2xy( const std::vector<double> &lat, const std::vector<double> &lon,
		std::vector<double> &x, std::vector<double> &y, double reflat ) {

	if (lat.size() != lon.size()) {
		throw std::runtime_error( "Mismatched input lat/lon vector lengths!" );
	}

	if (reflat < -90.0) {
		reflat = 0.5 * ( NCPA::max( lat ) - NCPA::min( lat ) );
	}
	double Re = NCPA::earthradius( reflat );
	x.clear();
	y.clear();
	x.resize( lat.size(), 0.0 );
	y.resize( lat.size(), 0.0 );

	for (size_t i = 0; i < lat.size(); i++) {
		x[ i ] = Re * NCPA::deg2rad( lon[ i ] ) * std::cos( NCPA::deg2rad( reflat ) );
		y[ i ] = Re * NCPA::deg2rad( lat[ i ] );
	}
}