#include "NCPACommon.h"
#include "NCPAAtmosphere.h"

#include <set>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <vector>


NCPA::Atmosphere3D::~Atmosphere3D() { }

void NCPA::Atmosphere3D::convert_1D_to_3D( const double *in, double ***out, size_t nz ) {
	out = NCPA::matrix3d<double>( 1, 1, nz );
	std::memcpy( out[0][0], in, nz*sizeof(double) );
}

// second derivative
double NCPA::Atmosphere3D::get_derivative( double x, double y,
		const std::string &key, deriv_t direction1,
		deriv_t direction2 ) const {
	deriv_t dvec[ 2 ] = { direction1, direction2 };
	return get_derivative( x, y, key, 2, dvec );
}

// second derivative
double NCPA::Atmosphere3D::get_derivative( double x, double y, double z,
		const std::string &key, deriv_t direction1,
		deriv_t direction2 ) const {
	deriv_t dvec[ 2 ] = { direction1, direction2 };
	return get_derivative( x, y, z, key, 2, dvec );
}

// first derivative
double NCPA::Atmosphere3D::get_derivative( double x, double y,
			const std::string &key, deriv_t direction ) const {
	return get_derivative( x, y, key, 1, &direction );
}

// first derivative
double NCPA::Atmosphere3D::get_derivative( double x, double y, double z,
				const std::string &key, deriv_t direction ) const {
	return get_derivative( x, y, z, key, 1, &direction );
}
