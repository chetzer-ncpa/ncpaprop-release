#include "AtmosphericProperty3D.h"
#include "LANLInterpolation.h"
#include "units.h"
#include "util.h"
#include <cstring>

// base virtual destructor
NCPA::AtmosphericProperty3D::~AtmosphericProperty3D() { }

NCPA::units_t NCPA::AtmosphericProperty3D::get_range_units() const { return r_units_; }

NCPA::units_t NCPA::AtmosphericProperty3D::get_property_units() const { return f_units_; }

NCPA::units_t NCPA::AtmosphericProperty3D::get_altitude_units() const { return z_units_; }


NCPA::VectorAtmosphericProperty3D::VectorAtmosphericProperty3D(
			const std::string &key,
			size_t nx, double *xvals,
			size_t ny, double *yvals,
			NCPA::units_t range_units,
			NCPA::AtmosphericProperty1D ***prop_mat,
			size_t nz, double *zvals, NCPA::units_t altitude_units ) {

	key_ = key;
	LANL::prep( spline_, nx, ny, nz );
	NCPA::AtmosphericProperty1D *temp_prop;

	std::cout << "Creating key " << key << std::endl;

	// populate index variables
	size_t i, j, k;
	for (i = 0; i < nx; i++) {
		spline_.x_vals[ i ] = xvals[ i ];
	}
	for (j = 0; j < ny; j++) {
		spline_.y_vals[ j ] = yvals[ j ];
	}
	for (k = 0; k < nz; k++) {
		spline_.z_vals[ k ] = zvals[ k ];
	}


	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			temp_prop = prop_mat[ i ][ j ];
			for (k = 0; k < nz; k++) {
				spline_.f_vals[ i ][ j ][ k ] = temp_prop->get( zvals[ k ] );
			}
		}
	}

	// Create the internal splines
	LANL::set( spline_ );
	r_units_ = range_units;
	z_units_ = altitude_units;
	f_units_ = prop_mat[ 0 ][ 0 ]->get_units();
}

NCPA::VectorAtmosphericProperty3D::VectorAtmosphericProperty3D(
	const std::string &key, size_t nx, double *xvals, size_t ny, double *yvals,
	size_t nz, double *zvals, double ***prop_mat, NCPA::units_t range_units,
	NCPA::units_t altitude_units, NCPA::units_t property_units ) {

	key_ = key;
	LANL::prep( spline_, nx, ny, nz );

	// populate index variables
	size_t i, j, k;
	for (i = 0; i < nx; i++) {
		spline_.x_vals[ i ] = xvals[ i ];
	}
	for (j = 0; j < ny; j++) {
		spline_.y_vals[ j ] = yvals[ j ];
	}
	for (k = 0; k < nz; k++) {
		spline_.z_vals[ k ] = zvals[ k ];
	}


	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				spline_.f_vals[ i ][ j ][ k ] = prop_mat[ i ][ j ][ k ];
			}
		}
	}

	// Create the internal splines
	LANL::set( spline_ );
	r_units_ = range_units;
	z_units_ = altitude_units;
	f_units_ = property_units;
}

// copy constructor
NCPA::VectorAtmosphericProperty3D::VectorAtmosphericProperty3D(
	const NCPA::VectorAtmosphericProperty3D &prop ) {

	int i, j, k;

	this->z_units_ = prop.z_units_;
	this->r_units_ = prop.r_units_;
	this->f_units_ = prop.f_units_;
	this->key_ = prop.key_;

	LANL::prep( this->spline_, prop.spline_.length_x, prop.spline_.length_y,
		prop.spline_.length_z );
	for (i = 0; i < prop.spline_.length_x; i++) {
		this->spline_.x_vals[ i ] = prop.spline_.x_vals[ i ];
	}
	for (j = 0; j < prop.spline_.length_y; j++) {
		this->spline_.y_vals[ j ] = prop.spline_.y_vals[ j ];
	}
	for (k = 0; k < prop.spline_.length_z; k++) {
		this->spline_.z_vals[ k ] = prop.spline_.z_vals[ k ];
	}
	for (i = 0; i < this->spline_.length_x; i++) {
		for (j = 0; j < this->spline_.length_y; j++) {
			for (k = 0; k < this->spline_.length_z; k++) {
				this->spline_.f_vals[ i ][ j ][ k ] = prop.spline_.f_vals[ i ][ j ][ k ];
			}
		}
	}

	LANL::set( this->spline_ );
}

NCPA::VectorAtmosphericProperty3D::~VectorAtmosphericProperty3D() {
	LANL::clear( spline_ );
}

double NCPA::VectorAtmosphericProperty3D::get( double x, double y, double z ) {
	return LANL::eval_f( x, y, z, spline_ );
}

double NCPA::VectorAtmosphericProperty3D::get( double x, double y ) {
	return 0.0;
}

void NCPA::VectorAtmosphericProperty3D::convert_range_units( NCPA::units_t new_units ) {
	// get original data vectors from spline
	double *x, *y, *z, ***f;
	size_t i, j;
	size_t nx = spline_.length_x, ny = spline_.length_y, nz = spline_.length_z;
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
	z = new double[ nz ];
	std::memcpy( z, spline_.z_vals, nz*sizeof(double) );
	f = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			std::memcpy( f[ i ][ j ], spline_.f_vals[ i ][ j ], nz*sizeof(double) );
		}
	}
	LANL::clear( spline_ );

	// convert
	NCPA::units_t old_units = r_units_;
	NCPA::Units::convert( x, nx, old_units, new_units, x );
	NCPA::Units::convert( y, ny, old_units, new_units, y );

	// regenerate spline
	LANL::prep( spline_, nx, ny, nz );
	std::memcpy( spline_.x_vals, x, nx*sizeof(double) );
	std::memcpy( spline_.y_vals, y, ny*sizeof(double) );
	std::memcpy( spline_.z_vals, z, nz*sizeof(double) );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			std::memcpy( spline_.f_vals[ i ][ j ], f[ i ][ j ], nz*sizeof(double) );
		}
	}
	LANL::set( spline_ );
	r_units_ = new_units;

	NCPA::free_matrix3d<double>( f, nx, ny, nz );
	delete [] x;
	delete [] y;
	delete [] z;
}

void NCPA::VectorAtmosphericProperty3D::convert_altitude_units( NCPA::units_t new_units ) {
	// get original data vectors from spline
	double *x, *y, *z, ***f;
	size_t i, j;
	size_t nx = spline_.length_x, ny = spline_.length_y, nz = spline_.length_z;
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
	z = new double[ nz ];
	std::memcpy( z, spline_.z_vals, nz*sizeof(double) );
	f = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			std::memcpy( f[ i ][ j ], spline_.f_vals[ i ][ j ], nz*sizeof(double) );
		}
	}
	LANL::clear( spline_ );

	// convert
	NCPA::units_t old_units = z_units_;
	NCPA::Units::convert( z, nz, old_units, new_units, z );

	// regenerate spline
	LANL::prep( spline_, nx, ny, nz );
	std::memcpy( spline_.x_vals, x, nx*sizeof(double) );
	std::memcpy( spline_.y_vals, y, ny*sizeof(double) );
	std::memcpy( spline_.z_vals, z, nz*sizeof(double) );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			std::memcpy( spline_.f_vals[ i ][ j ], f[ i ][ j ], nz*sizeof(double) );
		}
	}
	LANL::set( spline_ );
	z_units_ = new_units;

	NCPA::free_matrix3d<double>( f, nx, ny, nz );
	delete [] x;
	delete [] y;
	delete [] z;
}

void NCPA::VectorAtmosphericProperty3D::convert_property_units( NCPA::units_t new_units ) {
	// get original data vectors from spline
	double *x, *y, *z, ***f;
	size_t i, j;
	size_t nx = spline_.length_x, ny = spline_.length_y, nz = spline_.length_z;
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
	z = new double[ nz ];
	std::memcpy( z, spline_.z_vals, nz*sizeof(double) );
	f = NCPA::matrix3d<double>( nx, ny, nz );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			std::memcpy( f[ i ][ j ], spline_.f_vals[ i ][ j ], nz*sizeof(double) );
		}
	}
	LANL::clear( spline_ );

	// convert and regenerate spline
	NCPA::units_t old_units = f_units_;
	LANL::prep( spline_, nx, ny, nz );
	std::memcpy( spline_.x_vals, x, nx*sizeof(double) );
	std::memcpy( spline_.y_vals, y, ny*sizeof(double) );
	std::memcpy( spline_.z_vals, z, nz*sizeof(double) );
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			NCPA::Units::convert( f[ i ][ j ], nz, old_units, new_units, f[ i ][ j ] );
			std::memcpy( spline_.f_vals[ i ][ j ], f[ i ][ j ], nz*sizeof(double) );
		}
	}
	LANL::set( spline_ );
	f_units_ = new_units;

	NCPA::free_matrix3d<double>( f, nx, ny, nz );
	delete [] x;
	delete [] y;
	delete [] z;
}

double NCPA::VectorAtmosphericProperty3D::x_max() const {
	return spline_.x_vals[ spline_.length_x - 1 ];
}

double NCPA::VectorAtmosphericProperty3D::x_min() const {
	return spline_.x_vals[ 0 ];
}

double NCPA::VectorAtmosphericProperty3D::y_max() const {
	return spline_.y_vals[ spline_.length_y - 1 ];
}

double NCPA::VectorAtmosphericProperty3D::y_min() const {
	return spline_.y_vals[ 0 ];
}

double NCPA::VectorAtmosphericProperty3D::z_max() const {
	return spline_.z_vals[ spline_.length_z - 1 ];
}

double NCPA::VectorAtmosphericProperty3D::z_min() const {
	return spline_.z_vals[ 0 ];
}

size_t NCPA::VectorAtmosphericProperty3D::x_size() const {
	return (size_t)(spline_.length_x);
}

size_t NCPA::VectorAtmosphericProperty3D::y_size() const {
	return (size_t)(spline_.length_y);
}

size_t NCPA::VectorAtmosphericProperty3D::z_size() const {
	return (size_t)(spline_.length_z);
}

void NCPA::VectorAtmosphericProperty3D::x_vector( size_t &nx, double *&x ) const {
	nx = (size_t)(spline_.length_x);
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
}

void NCPA::VectorAtmosphericProperty3D::y_vector( size_t &ny, double *&y ) const {
	ny = (size_t)(spline_.length_y);
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
}

void NCPA::VectorAtmosphericProperty3D::z_vector( size_t &nz, double *&z ) const {
	nz = (size_t)(spline_.length_z);
	z = new double[ nz ];
	std::memcpy( z, spline_.z_vals, nz*sizeof(double) );
}

double NCPA::VectorAtmosphericProperty3D::get_derivative( double x, double y,
			deriv_t direction ) {
	throw std::runtime_error( "Scalar derivative requested for vector quantity "
		+ key_ + " (how did you manage that?)" );
}

double NCPA::VectorAtmosphericProperty3D::get_derivative( double x, double y,
			size_t order, deriv_t *directions ) {
	throw std::runtime_error( "Scalar derivative requested for vector quantity "
		+ key_ + " (how did you manage that?)" );
}

double NCPA::VectorAtmosphericProperty3D::get_derivative( double x, double y, double z,
			deriv_t direction ) {
	return get_derivative( x, y, z, 1, &direction );
}

double NCPA::VectorAtmosphericProperty3D::get_derivative( double x, double y, double z,
			size_t order, deriv_t *directions ) {
	switch (order) {
		case 1:
			return LANL::eval_df( x, y, z, (int)(directions[0]), spline_ );
		case 2:
			return LANL::eval_ddf( x, y, z, (int)(directions[0]), (int)(directions[1]),
				spline_ );
	}
	throw std::runtime_error( "Maximum second order derivative available" );
}

void NCPA::VectorAtmosphericProperty3D::as_matrix( double ***&data,
	size_t &nx, size_t &ny, size_t &nz ) const {

	int i, j, k;
	data = NCPA::matrix3d<double>( spline_.length_x, spline_.length_y, spline_.length_z );
	for (i = 0; i < spline_.length_x; i++) {
		for (j = 0; j < spline_.length_y; j++) {
			for (k = 0; k < spline_.length_z; k++) {
				data[ i ][ j ][ k ] = spline_.f_vals[ i ][ j ][ k ];
			}
		}
	}
	nx = spline_.length_x;
	ny = spline_.length_y;
	nz = spline_.length_z;
}

void NCPA::VectorAtmosphericProperty3D::free_matrix( double ***data ) const {
	NCPA::free_matrix3d<double>( data, spline_.length_x, spline_.length_y, spline_.length_z );
}



//////////////////////////////////////////////////////////////////////////////////


NCPA::ScalarAtmosphericProperty3D::ScalarAtmosphericProperty3D(
			const std::string &key,
			size_t nx, double *xvals,
			size_t ny, double *yvals,
			NCPA::units_t range_units,
			NCPA::ScalarWithUnits ***prop_mat ) {

	key_ = key;
	LANL::prep( spline_, nx, ny );
	NCPA::ScalarWithUnits *temp_prop;

	// populate index variables
	size_t i, j;
	for (i = 0; i < nx; i++) {
		spline_.x_vals[ i ] = xvals[ i ];
	}
	for (j = 0; j < ny; j++) {
		spline_.y_vals[ j ] = yvals[ j ];
	}


	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			temp_prop = prop_mat[ i ][ j ];
			spline_.f_vals[ i ][ j ] = prop_mat[ i ][ j ]->get();
		}
	}

	// Create the internal splines
	LANL::set( spline_ );
	r_units_ = range_units;
	f_units_ = prop_mat[ 0 ][ 0 ]->get_units();
	z_units_ = NCPA::UNITS_NONE;
}


NCPA::ScalarAtmosphericProperty3D::ScalarAtmosphericProperty3D(
	const std::string &key, size_t nx, double *xvals, size_t ny, double *yvals,
	double **prop_mat, NCPA::units_t range_units, NCPA::units_t property_units ) {

	key_ = key;
	LANL::prep( spline_, nx, ny );

	// populate index variables
	size_t i, j;
	for (i = 0; i < nx; i++) {
		spline_.x_vals[ i ] = xvals[ i ];
	}
	for (j = 0; j < ny; j++) {
		spline_.y_vals[ j ] = yvals[ j ];
	}


	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			spline_.f_vals[ i ][ j ] = prop_mat[ i ][ j ];
		}
	}

	// Create the internal splines
	LANL::set( spline_ );
	r_units_ = range_units;
	f_units_ = property_units;
	z_units_ = NCPA::UNITS_NONE;
}


NCPA::ScalarAtmosphericProperty3D::ScalarAtmosphericProperty3D(
	const NCPA::ScalarAtmosphericProperty3D &prop ) {

	int i, j;

	this->r_units_ = prop.r_units_;
	this->f_units_ = prop.f_units_;
	this->z_units_ = prop.z_units_;
	this->key_ = prop.key_;

	LANL::prep( this->spline_, prop.spline_.length_x, prop.spline_.length_y );
	for (i = 0; i < prop.spline_.length_x; i++) {
		this->spline_.x_vals[ i ] = prop.spline_.x_vals[ i ];
	}
	for (j = 0; j < prop.spline_.length_y; j++) {
		this->spline_.y_vals[ j ] = prop.spline_.y_vals[ j ];
	}
	for (i = 0; i < this->spline_.length_x; i++) {
		for (j = 0; j < this->spline_.length_y; j++) {
			this->spline_.f_vals[ i ][ j ] = prop.spline_.f_vals[ i ][ j ];
		}
	}
	LANL::set( this->spline_ );
}

NCPA::ScalarAtmosphericProperty3D::~ScalarAtmosphericProperty3D() {
	LANL::clear( spline_ );
}

double NCPA::ScalarAtmosphericProperty3D::get( double x, double y, double z ) {
	return LANL::eval_f( x, y, spline_ );
}

double NCPA::ScalarAtmosphericProperty3D::get( double x, double y ) {
	return LANL::eval_f( x, y, spline_ );
}

void NCPA::ScalarAtmosphericProperty3D::convert_range_units( NCPA::units_t new_units ) {
	// get original data vectors from spline
	double *x, *y, **f;
	size_t i;
	size_t nx = spline_.length_x, ny = spline_.length_y;
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
	f = NCPA::matrix<double>( nx, ny );
	for (i = 0; i < nx; i++) {
		std::memcpy( f[ i ], spline_.f_vals[ i ], ny*sizeof(double) );
	}
	LANL::clear( spline_ );

	// convert
	NCPA::units_t old_units = r_units_;
	NCPA::Units::convert( x, nx, old_units, new_units, x );
	NCPA::Units::convert( y, ny, old_units, new_units, y );

	// regenerate spline
	LANL::prep( spline_, nx, ny );
	std::memcpy( spline_.x_vals, x, nx*sizeof(double) );
	std::memcpy( spline_.y_vals, y, nx*sizeof(double) );
	for (i = 0; i < nx; i++) {
		std::memcpy( spline_.f_vals[ i ], f[ i ], ny*sizeof(double) );
	}
	LANL::set( spline_ );
	r_units_ = new_units;

	NCPA::free_matrix<double>( f, nx, ny );
	delete [] x;
	delete [] y;
}

void NCPA::ScalarAtmosphericProperty3D::convert_property_units( NCPA::units_t new_units ) {
	// get original data vectors from spline
	double *x, *y, **f;
	size_t i;
	size_t nx = spline_.length_x, ny = spline_.length_y;
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
	f = NCPA::matrix<double>( nx, ny );
	for (i = 0; i < nx; i++) {
		std::memcpy( f[ i ], spline_.f_vals[ i ], ny*sizeof(double) );
	}
	LANL::clear( spline_ );

	// convert and regenerate spline
	NCPA::units_t old_units = r_units_;
	LANL::prep( spline_, nx, ny );
	std::memcpy( spline_.x_vals, x, nx*sizeof(double) );
	std::memcpy( spline_.y_vals, y, nx*sizeof(double) );
	for (i = 0; i < nx; i++) {
		NCPA::Units::convert( f[ i ], ny, old_units, new_units, f[ i ] );
		std::memcpy( spline_.f_vals[ i ], f[ i ], ny*sizeof(double) );
	}
	LANL::set( spline_ );
	f_units_ = new_units;

	NCPA::free_matrix<double>( f, nx, ny );
	delete [] x;
	delete [] y;
}

void NCPA::ScalarAtmosphericProperty3D::convert_altitude_units(
	NCPA::units_t new_units ) {}

double NCPA::ScalarAtmosphericProperty3D::x_max() const {
	return spline_.x_vals[ spline_.length_x - 1 ];
}

double NCPA::ScalarAtmosphericProperty3D::x_min() const {
	return spline_.x_vals[ 0 ];
}

double NCPA::ScalarAtmosphericProperty3D::y_max() const {
	return spline_.y_vals[ spline_.length_y - 1 ];
}

double NCPA::ScalarAtmosphericProperty3D::y_min() const {
	return spline_.y_vals[ 0 ];
}

double NCPA::ScalarAtmosphericProperty3D::z_max() const {
	return 0.0;
}

double NCPA::ScalarAtmosphericProperty3D::z_min() const {
	return 0.0;
}

size_t NCPA::ScalarAtmosphericProperty3D::x_size() const {
	return (size_t)(spline_.length_x);
}

size_t NCPA::ScalarAtmosphericProperty3D::y_size() const {
	return (size_t)(spline_.length_y);
}

size_t NCPA::ScalarAtmosphericProperty3D::z_size() const {
	return 0;
}

void NCPA::ScalarAtmosphericProperty3D::x_vector( size_t &nx, double *&x ) const {
	nx = (size_t)(spline_.length_x);
	x = new double[ nx ];
	std::memcpy( x, spline_.x_vals, nx*sizeof(double) );
}

void NCPA::ScalarAtmosphericProperty3D::y_vector( size_t &ny, double *&y ) const {
	ny = (size_t)(spline_.length_y);
	y = new double[ ny ];
	std::memcpy( y, spline_.y_vals, ny*sizeof(double) );
}

void NCPA::ScalarAtmosphericProperty3D::z_vector( size_t &nz, double *&z ) const {
	nz = 0;
	z = NULL;
}

double NCPA::ScalarAtmosphericProperty3D::get_derivative( double x, double y, double z,
			size_t order, deriv_t *directions ) {
	throw std::runtime_error( "Vector derivative requested for scalar quantity "
		+ key_ + " (how did you manage that?)" );
}

double NCPA::ScalarAtmosphericProperty3D::get_derivative( double x, double y, double z,
			deriv_t direction ) {
	throw std::runtime_error( "Vector derivative requested for scalar quantity "
		+ key_ + " (how did you manage that?)" );
}

double NCPA::ScalarAtmosphericProperty3D::get_derivative( double x, double y,
			size_t order, deriv_t *directions ) {
	switch (order) {
		case 1:
			return LANL::eval_df( x, y, (int)(directions[0]), spline_ );
		case 2:
			return LANL::eval_ddf( x, y, (int)(directions[0]), (int)(directions[1]),
				spline_ );
	}
	throw std::runtime_error( "Maximum second order derivative available" );
}

double NCPA::ScalarAtmosphericProperty3D::get_derivative( double x, double y,
			deriv_t direction ) {
	return get_derivative( x, y, 1, &direction );
}

void NCPA::ScalarAtmosphericProperty3D::as_matrix( double ***&data,
	size_t &nx, size_t &ny, size_t &nz ) const {

	int i, j;
	data = NCPA::matrix3d<double>( spline_.length_x, spline_.length_y, 1 );
	for (i = 0; i < spline_.length_x; i++) {
		for (j = 0; j < spline_.length_y; j++) {
			data[ i ][ j ][ 0 ] = spline_.f_vals[ i ][ j ];
		}
	}
	nx = spline_.length_x;
	ny = spline_.length_y;
	nz = 1;
}

void NCPA::ScalarAtmosphericProperty3D::free_matrix( double ***data ) const {
	NCPA::free_matrix3d<double>( data, spline_.length_x, spline_.length_y, 1 );
}