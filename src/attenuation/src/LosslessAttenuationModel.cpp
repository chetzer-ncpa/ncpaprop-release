#include "AttenuationModel.h"
#include "LosslessAttenuationModel.h"
#include <algorithm>
#include <vector>

NCPA::LosslessAttenuationModel::LosslessAttenuationModel()
		: NCPA::AttenuationModel() {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
}

NCPA::LosslessAttenuationModel::LosslessAttenuationModel(
		const NCPA::LosslessAttenuationModel &other )
		: NCPA::AttenuationModel( other ) {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
}

NCPA::LosslessAttenuationModel::LosslessAttenuationModel(
		NCPA::LosslessAttenuationModel &&other )
		: NCPA::AttenuationModel() {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
	::swap( *this, other );
}

void swap( NCPA::LosslessAttenuationModel &a, NCPA::LosslessAttenuationModel &b ) {
	using std::swap;
	::swap( static_cast<NCPA::AttenuationModel&>(a),
			static_cast<NCPA::AttenuationModel&>(b) );
}

NCPA::LosslessAttenuationModel& NCPA::LosslessAttenuationModel::operator=(
		NCPA::LosslessAttenuationModel other ) {
	::swap(*this,other);
	return *this;
}

NCPA::AttenuationModel* NCPA::LosslessAttenuationModel::clone() const {
	return static_cast<NCPA::AttenuationModel *>( new NCPA::LosslessAttenuationModel( *this ) );
}

NCPA::LosslessAttenuationModel::~LosslessAttenuationModel() {}

double NCPA::LosslessAttenuationModel::attenuation(double f, double z, NCPA::units_t z_units) {
	return 0.0;
}

void NCPA::LosslessAttenuationModel::attenuation(size_t n, double f, double *z, double *alpha,
		NCPA::units_t z_units) {
	std::fill(alpha, alpha+n, 0.0);
}

NCPA::units_t NCPA::LosslessAttenuationModel::get_calculation_units( attenuation_parameter_t param ) const {
	return NCPA::units_t::NONE;
}

const std::vector<NCPA::attenuation_parameter_t>
NCPA::LosslessAttenuationModel::get_required_parameters() const {
	return std::vector<NCPA::attenuation_parameter_t>();
}

NCPA::attenuation_model_t NCPA::LosslessAttenuationModel::type() const {
	return NCPA::attenuation_model_t::LOSSLESS;
}

NCPA::AttenuationModel* NCPA::LosslessAttenuationModel::finalize() {
	this->finalized_ = true;
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::LosslessAttenuationModel::set_parameter(attenuation_parameter_t param,
		const NCPA::VectorWithUnits &p) {
	return static_cast<NCPA::AttenuationModel *>( this );
}
