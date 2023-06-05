#ifndef NCPA__SUTHERLANDBASSATTENUATIONMODEL_H_INCLUDED
#define NCPA__SUTHERLANDBASSATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"
#include <vector>

namespace NCPA {

	class SutherlandBassAttenuationModel : public AttenuationModel {
	public:
		SutherlandBassAttenuationModel() {}
		virtual ~SutherlandBassAttenuationModel() {}

		virtual double attenuation(double f, double z, NCPA::units_t z_units = UNITS_NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units = UNITS_NONE);
		virtual units_t get_calculation_units( attenuation_model_parameter_label_t param ) const;
		virtual const std::vector<attenuation_model_parameter_label_t> get_required_parameters() const;

	protected:

		const double mu_o  		= 18.192E-6;    		// Reference viscosity [kg/(m*s)]
		const double T_o   		= 293.15;         		// Reference temperature [K]
		const double P_o   		= 101325;        		// Reference pressure [Pa]
		const double S     		= 117.0;	     		// Sutherland constant [K]
		const double gamma		= 1.4;

		// heat capacity|volume for O2, N2, CO2, and O3
		const double Cv_R[4] = { 5.0/2.0, 5.0/2.0, 3.0, 3.0 };
		// heat capacity|pressure for O2, N2, CO2, and O3
		const double Cp_R[4] = { 7.0/2.0, 7.0/2.0, 4.0, 4.0 };
		// characteristic temperature for O2, N2, CO2, and O3
		const double theta[4] = { 2239.1, 3352.0, 915.0, 1037.0 };

		const std::vector<attenuation_model_parameter_label_t> required_parameters = {
				ATTENUATION_MODEL_PARAMETER_ALTITUDE,
				ATTENUATION_MODEL_PARAMETER_TEMPERATURE,
				ATTENUATION_MODEL_PARAMETER_DENSITY,
				ATTENUATION_MODEL_PARAMETER_PRESSURE
		};


	}; // class SutherlandBassAttenuationModel

} // namespace

#endif
