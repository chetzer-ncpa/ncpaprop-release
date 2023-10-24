#include "AtmosphericProperty1D.h"
#include "NCPAUnits.h"
#include "NCPACommon.h"
#include "NCPAInterpolation.h"
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
//
//#include "gsl/gsl_errno.h"
//#include "gsl/gsl_spline.h"
//#include "gsl/gsl_version.h"

NCPA::AtmosphericProperty1D::AtmosphericProperty1D() : NCPA::vectorpair_t(),
	interp_type_{NCPAPROP_DEFAULT_ATMOSPHERIC_PROPERTY_INTERPOLATION_TYPE} {

	interp_ = NCPA::Interpolator1D::build( interp_type_ );
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D(
		NCPA::VectorWithUnits zvector,
		NCPA::VectorWithUnits propvector,
		NCPA::interpolator1d_t interptype ) : NCPA::vectorpair_t(zvector,propvector),
				interp_type_{interptype} {
	interp_ = NCPA::Interpolator1D::build( interp_type_ );
	build_splines_();
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units,
			NCPA::interpolator1d_t interptype ) : NCPA::vectorpair_t(),
					interp_type_{interptype} {
	this->first.set( n_points, altitude_points, altitude_units );
	this->second.set( n_points, property_values, property_units );
	interp_ = NCPA::Interpolator1D::build( interp_type_ );
	build_splines_();
}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points,
		double *altitude_points, const std::string &altitude_units,
		double *property_values, const std::string &property_units,
		NCPA::interpolator1d_t interptype )
			: NCPA::AtmosphericProperty1D(
				n_points, altitude_points, NCPA::Units::fromString(altitude_units),
				property_values, NCPA::Units::fromString(property_units)) {}

NCPA::AtmosphericProperty1D::AtmosphericProperty1D( const AtmosphericProperty1D &source )
		: NCPA::vectorpair_t(source) {
	this->first = source.first; // test that copy constructor is properly called
	this->second = source.second;
	this->interp_type_ = source.interp_type_;
	this->interp_ = NCPA::Interpolator1D::build( source.interp_type_ );
	reset_splines();
}

NCPA::AtmosphericProperty1D::~AtmosphericProperty1D() {
	delete_splines_();
}

void swap( NCPA::AtmosphericProperty1D &a, NCPA::AtmosphericProperty1D &b ) noexcept {
	using std::swap;
	swap(static_cast<NCPA::vectorpair_t&>(a), static_cast<NCPA::vectorpair_t&>(b));
	swap(a.interp_,b.interp_);
	swap(a.interp_type_,b.interp_type_);
}

NCPA::AtmosphericProperty1D &NCPA::AtmosphericProperty1D::operator=( NCPA::AtmosphericProperty1D other ) {
	::swap(*this,other);
	return *this;
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_altitude_units() const {
	NCPA::VectorWithUnits v(this->first);
	return v.get_units();
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_units() const {
	NCPA::VectorWithUnits v(this->second);
	return v.get_units();
}

size_t NCPA::AtmosphericProperty1D::size() const {
	return this->first.size();
}

void NCPA::AtmosphericProperty1D::convert_altitude_units( NCPA::units_t new_units ) {
	if (new_units != this->first.get_units()) {
		this->first.convert_units( new_units );
		this->reset_splines();
	}
}

void NCPA::AtmosphericProperty1D::convert_altitude_units( const std::string &new_units ) {
	this->convert_altitude_units( NCPA::Units::fromString(new_units) );
}

void NCPA::AtmosphericProperty1D::convert_units( NCPA::units_t new_units ) {
	this->second.convert_units( new_units );
	this->reset_splines();
}

void NCPA::AtmosphericProperty1D::convert_units( const std::string &new_units ) {
	this->convert_units( NCPA::Units::fromString(new_units) );
}

void NCPA::AtmosphericProperty1D::build_splines_() {
	if (interp_ == nullptr) {
		interp_ = NCPA::Interpolator1D::build(this->interp_type_);
	}

	double *zbuffer, *xbuffer;
	this->first.as_array(zbuffer);
	this->second.as_array(xbuffer);
	interp_->init()->allocate(this->first.size())->set(
			this->first.size(), zbuffer, xbuffer )->ready();
	delete [] zbuffer;
	delete [] xbuffer;
}

void NCPA::AtmosphericProperty1D::delete_splines_() {
	if (interp_ != nullptr) {
		interp_->free();
		delete interp_;
		interp_ = nullptr;
	}
}

void NCPA::AtmosphericProperty1D::reset_splines() {
	this->delete_splines_();
	this->build_splines_();
}

void NCPA::AtmosphericProperty1D::get_altitude_vector(
		double *buffer, units_t &buffer_units ) {
	this->first.get_values( buffer );
	buffer_units = this->first.get_units();
}

void NCPA::AtmosphericProperty1D::get_altitude_vector( double *buffer ) {
	this->first.get_values( buffer );
}

NCPA::VectorWithUnits NCPA::AtmosphericProperty1D::get_altitude_vector() {
	NCPA::VectorWithUnits v( this->first );
	return v;
}

void NCPA::AtmosphericProperty1D::get_altitude_vector_as(
		double *buffer, units_t buffer_units ) {
	NCPA::VectorWithUnits v( this->first );
	v.convert_units(buffer_units);
	v.get_values( buffer );
}

void NCPA::AtmosphericProperty1D::get_altitude_vector_as(
		double *buffer, const std::string &buffer_units ) {
	this->get_altitude_vector_as( buffer, NCPA::Units::fromString( buffer_units ) );
}


void NCPA::AtmosphericProperty1D::as_array( double *&buffer, units_t &buffer_units ) {
	this->second.as_array( buffer, buffer_units );
}

void NCPA::AtmosphericProperty1D::as_array( double *&buffer ) {
	this->second.as_array( buffer );
}

void NCPA::AtmosphericProperty1D::as_array( NCPA::ScalarWithUnits *&buffer ) {
	this->second.as_array( buffer );
}

NCPA::VectorWithUnits NCPA::AtmosphericProperty1D::get_vector() {
//	NCPA::VectorWithUnits v( this->second );
	return NCPA::VectorWithUnits( this->second );
}

NCPA::VectorWithUnits NCPA::AtmosphericProperty1D::get_vector_as( units_t buffer_units ) {
	NCPA::VectorWithUnits v( this->second );
	v.convert_units(buffer_units);
	return v;
}

NCPA::VectorWithUnits NCPA::AtmosphericProperty1D::get_vector_as( const std::string &buffer_units ) {
	return this->get_vector_as( NCPA::Units::fromString( buffer_units ) );
}

int NCPA::AtmosphericProperty1D::check_altitude_( double z_req ) const {
	if ( z_req < this->first.front() ) {
		return -1;
	} else if ( z_req > this->first.back() ) {
		return 1;
	} else {
		return 0;
	}
}

double NCPA::AtmosphericProperty1D::get( double z_req ) const {
	int check = check_altitude_( z_req );
	if ( check > 0 ) {
		return this->second.back().get();
	} else if ( check < 0 ) {
		return this->second.front().get();
	} else {
		return interp_->f( z_req );
//		return gsl_spline_eval( spline_, z_req, accel_ );
	}
}

double NCPA::AtmosphericProperty1D::get_as( double z_req, NCPA::units_t as_units ) const {
	return NCPA::Units::convert( this->get( z_req ), this->get_units(), as_units );
}

double NCPA::AtmosphericProperty1D::get_as( double z_req, const std::string &as_units ) const {
	return NCPA::Units::convert( this->get( z_req ), this->get_units(), as_units );
}

double NCPA::AtmosphericProperty1D::get_first_derivative( double z_req ) const {
	int check = check_altitude_( z_req );
	if (check == 0) {
		return interp_->df( z_req );
//		return gsl_spline_eval_deriv( spline_, z_req, accel_ );
	} else {
		return 0;
	}
}

double NCPA::AtmosphericProperty1D::get_first_derivative_as( double z_req, NCPA::units_t as_units ) const {
	return NCPA::Units::convert( this->get_first_derivative( z_req ), this->get_units(), as_units );
}

double NCPA::AtmosphericProperty1D::get_first_derivative_as( double z_req, const std::string &as_units ) const {
	return NCPA::Units::convert( this->get_first_derivative( z_req ), this->get_units(), as_units );
}

double NCPA::AtmosphericProperty1D::get_second_derivative( double z_req ) const {
	int check = check_altitude_( z_req );
	if (check == 0) {
		return interp_->df( 2, z_req );
//		return gsl_spline_eval_deriv2( spline_, z_req, accel_ );
	} else {
		return 0;
	}
}

double NCPA::AtmosphericProperty1D::get_second_derivative_as( double z_req, NCPA::units_t as_units ) const {
	return NCPA::Units::convert( this->get_second_derivative( z_req ), this->get_units(), as_units );
}

double NCPA::AtmosphericProperty1D::get_second_derivative_as( double z_req, const std::string &as_units ) const {
	return NCPA::Units::convert( this->get_second_derivative( z_req ), this->get_units(), as_units );
}

void NCPA::AtmosphericProperty1D::resample( double new_dz ) {

	double z0 = this->first.front().get();
	double z1 = this->first.back().get();
	size_t new_nz = (size_t)std::floor( ( z1 - z0 ) / new_dz ) + 1;

	double *new_z = new double[ new_nz ];
	double *new_prop = new double[ new_nz ];
	for (size_t i = 0; i < new_nz; i++) {
		new_z[ i ] = z0 + ((double)i) * new_dz;
		new_prop[ i ] = this->get( new_z[ i ] );
	}

	this->first.set( new_nz, new_z, this->get_altitude_units() );
	this->second.set( new_nz, new_prop, this->get_units() );
	this->reset_splines();
	delete [] new_z;
	delete [] new_prop;
}
