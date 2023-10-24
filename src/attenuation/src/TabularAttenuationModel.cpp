#include "AttenuationModel.h"
#include "TabularAttenuationModel.h"
#include "NCPAUnits.h"
#include "NCPAInterpolation.h"


NCPA::TabularAttenuationModel::TabularAttenuationModel(
		const NCPA::TabularAttenuationModel &other )
		: NCPA::AttenuationModel( other ) {}

NCPA::TabularAttenuationModel::TabularAttenuationModel(
		NCPA::TabularAttenuationModel &&other )
		: NCPA::AttenuationModel() {
	::swap( *this, other );
}

void swap( NCPA::TabularAttenuationModel &a, NCPA::TabularAttenuationModel &b ) {
	using std::swap;
	::swap( static_cast<NCPA::AttenuationModel&>(a),
			static_cast<NCPA::AttenuationModel&>(b) );
}

NCPA::TabularAttenuationModel& NCPA::TabularAttenuationModel::operator=(
		NCPA::TabularAttenuationModel other ) {
	::swap(*this,other);
	return *this;
}

NCPA::AttenuationModel* NCPA::TabularAttenuationModel::clone() const {
	return static_cast<NCPA::AttenuationModel *>( new NCPA::TabularAttenuationModel( *this ) );
}

NCPA::TabularAttenuationModel::~TabularAttenuationModel() {}

NCPA::AttenuationModel* NCPA::TabularAttenuationModel::set_parameter(attenuation_parameter_t param,
		const NCPA::VectorWithUnits &p) {
	if (z_vector_.empty()) {
		throw std::out_of_range( "AttenuationModel: No existing vertical scale set, so can't interpret new vertical parameter." );
	}
	if (z_vector_.size() != p.size()) {
		throw std::out_of_range( "AttenuationModel: Existing vertical scale has different number of points than new vertical parameter." );
	}
	this->set_parameter_(param, p);
	if (param == NCPA::attenuation_parameter_t::ATTENUATION) {
		native_output_units_ = this->vector_parameters_[param].get_units();
	}
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::TabularAttenuationModel::set_parameter(attenuation_parameter_t param,
		size_t nz, const double *p, NCPA::units_t units) {
//	if (z_vector_.empty()) {
//		throw std::out_of_range( "AttenuationModel: No existing vertical scale set, so can't interpret new vertical parameter." );
//	}
//	if (z_vector_.size() != nz) {
//		throw std::out_of_range( "AttenuationModel: Existing vertical scale has different number of points than new vertical parameter." );
//	}
//	if (param == NCPA::attenuation_parameter_t::ATTENUATION) {
//		native_output_units_ = units;
//	}
//	NCPA::VectorWithUnits pp( nz, p, units );
	return this->set_parameter( param, NCPA::VectorWithUnits( nz, p, units ) );
//	this->set_parameter_(param, pp);
//	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::TabularAttenuationModel::set_parameter(attenuation_parameter_t param,
		const NCPA::ScalarWithUnits &p) {
	this->set_parameter_(param, p);
	if (param == NCPA::attenuation_parameter_t::ATTENUATION) {
		native_output_units_ = this->scalar_parameters_[param].get_units();
	}
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::TabularAttenuationModel::set_parameter(attenuation_parameter_t param,
		double p, NCPA::units_t units) {
//	if (param == NCPA::attenuation_parameter_t::ATTENUATION) {
//		native_output_units_ = units;
//	}
//	NCPA::ScalarWithUnits s( p, units );
//	this->set_parameter_(param,s);
	return this->set_parameter( param, NCPA::ScalarWithUnits( p, units ) );
}

void NCPA::TabularAttenuationModel::attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units) {
	for (size_t i = 0; i < n; i++) {
		double zconv = NCPA::Units::convert( z[i], z_units, z_vector_.get_units() );
		alpha[i] = NCPA::Units::convert(
				interpolators_[NCPA::attenuation_parameter_t::ATTENUATION]->f(zconv),
				native_output_units_,
				output_units_ );
	}
}

double NCPA::TabularAttenuationModel::attenuation(double f, double z, NCPA::units_t z_units ) {
	double a;
	this->attenuation( 1, f, &z, &a, z_units );
	return a;
}

NCPA::units_t NCPA::TabularAttenuationModel::get_calculation_units(
		NCPA::attenuation_parameter_t param ) const {
	return NCPA::units_t::NONE;
}

const std::vector<NCPA::attenuation_parameter_t>
NCPA::TabularAttenuationModel::get_required_parameters() const {
	return std::vector<NCPA::attenuation_parameter_t>( required_parameters_ );
}

NCPA::attenuation_model_t NCPA::TabularAttenuationModel::type() const {
	return NCPA::attenuation_model_t::TABULAR;
}




