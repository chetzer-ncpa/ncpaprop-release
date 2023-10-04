#include "ISO9613_1AttenuationModel.h"
#include "NCPAUnits.h"
#include <cmath>
#include <vector>


NCPA::units_t NCPA::ISO9613_1AttenuationModel::get_calculation_units(
		attenuation_model_parameter_label_t param ) const {
	switch (param) {
	case ATTENUATION_MODEL_PARAMETER_ALTITUDE:		// default km
		return NCPA::UNITS_DISTANCE_KILOMETERS;
		break;
	case ATTENUATION_MODEL_PARAMETER_TEMPERATURE:	// default K
		return NCPA::UNITS_TEMPERATURE_KELVIN;
		break;
	case ATTENUATION_MODEL_PARAMETER_PRESSURE:		// default atm
		return NCPA::UNITS_PRESSURE_ATMOSPHERES;
		break;
	case ATTENUATION_MODEL_PARAMETER_HUMIDITY:		// default Pa
		return NCPA::UNITS_NONE;
		break;
	default:
		return NCPA::UNITS_NONE;
	}
}

const std::vector<NCPA::attenuation_model_parameter_label_t>
NCPA::ISO9613_1AttenuationModel::get_required_parameters() const {
	return required_parameters;
}

double NCPA::ISO9613_1AttenuationModel::attenuation(
		double f, double z, NCPA::units_t z_units) {
	double alpha;
	this->attenuation( 1, f, &z, &alpha, z_units );
	return alpha;
}

void NCPA::ISO9613_1AttenuationModel::attenuation(
		size_t n, double freq, double *zvec, double *alpha, NCPA::units_t z_units) {
	if (z_units == NCPA::UNITS_NONE) {
		z_units = this->get_calculation_units(ATTENUATION_MODEL_PARAMETER_ALTITUDE);
	}

	// get units right
	for (auto it = get_required_parameters().cbegin(); it != get_required_parameters().cend(); ++it) {
		this->parameter_map[*it]->convert_units( this->get_calculation_units( *it ) );
		this->parameter_map[*it]->convert_altitude_units( this->get_calculation_units(ATTENUATION_MODEL_PARAMETER_ALTITUDE) );
	}
//	z_vector_.convert_altitude_units( this->get_calculation_units(ATTENUATION_MODEL_PARAMETER_ALTITUDE) );

	for (size_t i = 0; i < n; i++) {
		double z = NCPA::Units::convert(zvec[i],z_units,
				this->get_calculation_units(ATTENUATION_MODEL_PARAMETER_ALTITUDE) );
		double P_z = this->parameter_map[ATTENUATION_MODEL_PARAMETER_PRESSURE]->get(z);
		double H_z = this->parameter_map[ATTENUATION_MODEL_PARAMETER_HUMIDITY]->get(z);
		double T_z = this->parameter_map[ATTENUATION_MODEL_PARAMETER_TEMPERATURE]->get(z);

		double F = freq / P_z;                // frequency per atm

		// calculate saturation pressure
		double Psat = std::pow( 10.0,
			(10.79586*(1-(To1/T_z)) - 5.02808*std::log10(T_z/To1)
				+ 1.50474e-4*(1.0 - std::pow(10.0, -8.29692*((T_z/To1)-1)))
				-  4.2873e-4*(1.0 - std::pow(10.0, -4.76955*((To1/T_z)-1)))
				- 2.2195983) );

		// absolute humidity
		double h = H_z * Psat / P_z;

		// scaled relaxation frequency for nitrogen
		double FrN = std::sqrt(To/T_z)
			* (9.0 + 280.0*h*std::exp(-4.17*(std::pow(To/T_z,1.0/3.0)-1.0)));

		// scaled relaxation frequency for oxygen
		double FrO = 24.0 + 4.04e4 * h * (h+0.02)/(h+0.391);

		// attenuation coefficient in nepers/m
		double term1 = P_z*F*F*(1.84e-11*std::sqrt(T_z/To));
		double term2 = std::pow(T_z/To,-2.5)
			* (1.275e-2 * std::exp(-2239.1/T_z) / (FrO + F*F/FrO));
		double term3 = 1.068e-1 * std::exp(-3352.0 / T_z) / (FrN + F*F/FrN);
		alpha[i] = term1 + term2 + term3;
	}
}
