#ifndef NCPA__ISO9613_1ATTENUATIONMODEL_H_INCLUDED
#define NCPA__ISO9613_1ATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"
#include <vector>

namespace NCPA {

	class ISO9613_1AttenuationModel : public AttenuationModel {
	public:
		ISO9613_1AttenuationModel() {}
		virtual ~ISO9613_1AttenuationModel() {}

		virtual double attenuation(double f, double z, NCPA::units_t z_units = UNITS_NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units = UNITS_NONE);
		virtual units_t get_calculation_units( attenuation_model_parameter_label_t param ) const;
		virtual const std::vector<attenuation_model_parameter_label_t> get_required_parameters() const;

	protected:

		const double To1 = 273.15;      // triple point
		const double To  = 293.15;      // reference Temperature

		const std::vector<attenuation_model_parameter_label_t> required_parameters = {
				ATTENUATION_MODEL_PARAMETER_TEMPERATURE,
				ATTENUATION_MODEL_PARAMETER_HUMIDITY,
				ATTENUATION_MODEL_PARAMETER_PRESSURE
		};


	}; // class ISO9613_1AttenuationModel

} // namespace

#endif
