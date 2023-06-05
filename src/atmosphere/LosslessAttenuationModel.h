#ifndef NCPA__ATTENUATION_LOSSLESSATTENUATIONMODEL_H_INCLUDED
#define NCPA__ATTENUATION_LOSSLESSATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"

namespace NCPA {
	class LosslessAttenuationModel : public AttenuationModel {
	public:
		LosslessAttenuationModel() {}
		virtual ~LosslessAttenuationModel() {}

		virtual double attenuation(double f, double z, NCPA::units_t z_units = UNITS_NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha, NCPA::units_t z_units = UNITS_NONE);
		virtual units_t get_calculation_units( attenuation_model_parameter_label_t param ) const;
		virtual const std::vector<attenuation_model_parameter_label_t> get_required_parameters() const;
	};
}

#endif
