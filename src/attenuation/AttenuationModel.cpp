#include <vector>

#include "AttenuationModel.h"
#include "NCPAUnits.h"

// Base class methods
NCPA::AttenuationModel::~AttenuationModel() {
	for (attenuation_parameter_map_t::iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
		delete it->second;
	}
}

void NCPA::AttenuationModel::set_vertical_scale(size_t nz, const double *z_km, NCPA::units_t units) {
	NCPA::VectorWithUnits *zvec = new NCPA::VectorWithUnits(nz, z_km, units );
	this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, zvec);
}

void NCPA::AttenuationModel::set_vertical_scale(size_t nz, const NCPA::ScalarWithUnits *z) {
	NCPA::VectorWithUnits *zvec = new NCPA::VectorWithUnits(nz, z);
	this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, zvec);
}

void NCPA::AttenuationModel::set_vertical_scale(const NCPA::VectorWithUnits *z) {
	NCPA::VectorWithUnits *zvec = new NCPA::VectorWithUnits(z);
	this->_set_parameter(ATTENUATION_MODEL_PARAMETER_ALTITUDE, zvec);
}

void NCPA::AttenuationModel::_set_parameter(attenuation_model_parameter_t param, const NCPA::VectorWithUnits *p) {
	attenuation_parameter_map_t::iterator it = parameter_map.find(param);
	if (it != parameter_map.end()) {
		delete it->second;
	}
	parameter_map[ param ] = p;
}
