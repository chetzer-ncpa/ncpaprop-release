#include <vector>
#include <stdexcept>

#include "AttenuationModel.h"
#include "AtmosphericProperty1D.h"
#include "NCPAUnits.h"

// Base class methods
NCPA::AttenuationModel::~AttenuationModel() {
	for (auto it = parameter_map.begin(); it != parameter_map.end(); ++it) {
		delete it->second;
	}
}

void NCPA::AttenuationModel::set_vertical_scale(const NCPA::VectorWithUnits &z) {
	z_vector_ = z;
	//NCPA::VectorWithUnits *zvec = new NCPA::VectorWithUnits(z);
	//this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, zvec);
}

void NCPA::AttenuationModel::set_vertical_scale(size_t nz, const double *p,
		NCPA::units_t units) {
	z_vector_ = NCPA::VectorWithUnits( nz, p, units );
//	this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, pp);
}

void NCPA::AttenuationModel::set_vertical_parameter(attenuation_model_parameter_label_t param,
		const NCPA::VectorWithUnits &p) {
	if (z_vector_.empty()) {
		throw std::out_of_range( "AttenuationModel: No existing vertical scale set, so can't interpret new vertical parameter." );
	}
	if (z_vector_.size() != p.size()) {
		throw std::out_of_range( "AttenuationModel: Existing vertical scale has different number of points than new vertical parameter." );
	}
	attenuation_model_parameter_t newparam = new NCPA::AtmosphericProperty1D( z_vector_, p );
	this->set_parameter_(param, newparam);
}

void NCPA::AttenuationModel::set_vertical_parameter(attenuation_model_parameter_label_t param,
		size_t nz, const double *p, NCPA::units_t units) {
	if (z_vector_.empty()) {
		throw std::out_of_range( "AttenuationModel: No existing vertical scale set, so can't interpret new vertical parameter." );
	}
	if (z_vector_.size() != nz) {
		throw std::out_of_range( "AttenuationModel: Existing vertical scale has different number of points than new vertical parameter." );
	}
	NCPA::VectorWithUnits pp( nz, p, units );
	attenuation_model_parameter_t newparam = new NCPA::AtmosphericProperty1D( z_vector_, pp );
	this->set_parameter_(param, newparam);
}

// this can have its own vertical scale because we will resample
void NCPA::AttenuationModel::set_vertical_parameter(attenuation_model_parameter_label_t param,
					const attenuation_model_parameter_t p) {
	this->set_parameter_( param, new NCPA::AtmosphericProperty1D( *p ) );
}

void NCPA::AttenuationModel::set_parameter_(
		attenuation_model_parameter_label_t param,
		attenuation_model_parameter_t p) {
	attenuation_parameter_map_t::iterator it = parameter_map.find(param);
	if (it != parameter_map.end()) {
		delete it->second;
	}
	parameter_map[ param ] = p;
}

