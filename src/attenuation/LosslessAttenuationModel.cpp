#include "LosslessAttenuationModel.h"
#include <algorithm>
#include <vector>

double NCPA::LosslessAttenuationModel::attenuation(double f, double z, NCPA::units_t z_units) {
	return 0.0;
}

void NCPA::LosslessAttenuationModel::attenuation(size_t n, double f, double *z, double *alpha,
		NCPA::units_t z_units) {
	std::fill(alpha, alpha+n, 0.0);
}

NCPA::units_t NCPA::LosslessAttenuationModel::get_calculation_units( attenuation_model_parameter_label_t param ) const {
	return NCPA::UNITS_NONE;
}

std::vector<NCPA::attenuation_model_parameter_label_t>
NCPA::LosslessAttenuationModel::get_required_parameters() const {
	return std::vector<NCPA::attenuation_model_parameter_label_t>();
}
