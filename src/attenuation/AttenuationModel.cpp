#include <vector>

#include "AttenuationModel.h"
#include "NCPAUnits.h"

// Base class methods
NCPA::AttenuationModel::~AttenuationModel() {
	for (attenuation_parameter_map_t::iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
		delete it->second;
	}
}

void NCPA::AttenuationModel::set_vertical_parameter(attenuation_model_parameter_label_t param,
		const NCPA::VectorWithUnits &p) {
	NCPA::VectorWithUnits *pp = new NCPA::VectorWithUnits( p );
	this->_set_parameter(param, pp);
}

void NCPA::AttenuationModel::set_vertical_parameter(attenuation_model_parameter_label_t param, size_t nz, const double *p,
		NCPA::units_t units) {
	NCPA::VectorWithUnits *pp = new NCPA::VectorWithUnits( nz, p, units );
	this->_set_parameter(param, pp);
}

void NCPA::AttenuationModel::set_vertical_scale(const NCPA::VectorWithUnits &z) {
	NCPA::VectorWithUnits *zvec = new NCPA::VectorWithUnits(z);
	this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, zvec);
}

void NCPA::AttenuationModel::set_vertical_scale(size_t nz, const double *p,
		NCPA::units_t units) {
	NCPA::VectorWithUnits *pp = new NCPA::VectorWithUnits( nz, p, units );
	this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, pp);
}

void NCPA::AttenuationModel::_set_parameter(attenuation_model_parameter_label_t param, NCPA::VectorWithUnits *p) {
	attenuation_parameter_map_t::iterator it = parameter_map.find(param);
	if (it != parameter_map.end()) {
		delete it->second;
	}
	parameter_map[ param ] = p;
}

