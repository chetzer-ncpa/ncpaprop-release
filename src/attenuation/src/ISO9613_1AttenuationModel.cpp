/**
 * Derived from:
 *
 * International Organization for Standardization.  1993.  Acoustics -
 *   Attenuation of sound during propagation outdoors - Part 1: Calculation
 *   of the absorption of sound by the atmosphere (ISO Standard No.
 *   9613-1:1993(E)).
 *
 */

#include "ISO9613_1AttenuationModel.h"
#include "NCPAUnits.h"
#include <cmath>
#include <vector>


NCPA::ISO9613_1AttenuationModel::ISO9613_1AttenuationModel()
		: NCPA::AttenuationModel() {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
}

NCPA::ISO9613_1AttenuationModel::ISO9613_1AttenuationModel(
		const NCPA::ISO9613_1AttenuationModel &other )
		: NCPA::AttenuationModel( other ) {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
}

NCPA::ISO9613_1AttenuationModel::ISO9613_1AttenuationModel(
		NCPA::ISO9613_1AttenuationModel &&other )
		: NCPA::AttenuationModel() {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
	::swap( *this, other );
}

void swap( NCPA::ISO9613_1AttenuationModel &a, NCPA::ISO9613_1AttenuationModel &b ) {
	using std::swap;
	::swap( static_cast<NCPA::AttenuationModel&>(a),
			static_cast<NCPA::AttenuationModel&>(b) );

}

NCPA::ISO9613_1AttenuationModel& NCPA::ISO9613_1AttenuationModel::operator=(
		NCPA::ISO9613_1AttenuationModel other ) {
	::swap(*this,other);
	return *this;
}

NCPA::AttenuationModel* NCPA::ISO9613_1AttenuationModel::clone() const {
	return static_cast<NCPA::AttenuationModel *>( new NCPA::ISO9613_1AttenuationModel( *this ) );
}

NCPA::ISO9613_1AttenuationModel::~ISO9613_1AttenuationModel() {}

NCPA::units_t NCPA::ISO9613_1AttenuationModel::get_calculation_units(
		attenuation_parameter_t param ) const {
	switch (param) {
	case NCPA::attenuation_parameter_t::ALTITUDE:		// default km
		return NCPA::units_t::DISTANCE_KILOMETERS;
		break;
	case NCPA::attenuation_parameter_t::TEMPERATURE:	// default K
		return NCPA::units_t::TEMPERATURE_KELVIN;
		break;
	case NCPA::attenuation_parameter_t::PRESSURE:		// default atm
		return NCPA::units_t::PRESSURE_ATMOSPHERES;
		break;
	case NCPA::attenuation_parameter_t::HUMIDITY:		// default n/a
		return NCPA::units_t::RATIO_PERCENT;
		break;
	default:
		return NCPA::units_t::NONE;
	}
}

const std::vector<NCPA::attenuation_parameter_t>
NCPA::ISO9613_1AttenuationModel::get_required_parameters() const {
	std::vector<attenuation_parameter_t> v(required_parameters_);
	return v;
}

double NCPA::ISO9613_1AttenuationModel::getval_(
		NCPA::attenuation_parameter_t param, double z ) {
	if (this->interpolators_.find(param) != this->interpolators_.end()) {
		return this->interpolators_[param]->f(z);
	} else if (this->scalar_parameters_.find(param) != this->scalar_parameters_.end()) {
		return this->scalar_parameters_[param].get();
	} else {
		std::ostringstream oss;
		oss << "No parameter set for " << NCPA::AttenuationModel::as_string(param);
		throw std::out_of_range(oss.str());
	}
}

double NCPA::ISO9613_1AttenuationModel::attenuation(
		double f, double z, NCPA::units_t z_units) {
	double alpha;
	this->attenuation( 1, f, &z, &alpha, z_units );
	return alpha;
}

void NCPA::ISO9613_1AttenuationModel::attenuation(
		size_t n, double freq, double *zvec, double *alpha, NCPA::units_t z_units) {
	if (!this->finalized_) {
		throw std::logic_error("Attenuation model not finalized!");
	}
	if (z_units == NCPA::units_t::NONE) {
		z_units = this->get_calculation_units(NCPA::attenuation_parameter_t::ALTITUDE);
	}

	for (size_t i = 0; i < n; i++) {
		double z = NCPA::Units::convert(zvec[i],z_units,
				this->get_calculation_units(NCPA::attenuation_parameter_t::ALTITUDE) );
		double P_z = this->getval_( NCPA::attenuation_parameter_t::PRESSURE, z );
		double H_z = this->getval_( NCPA::attenuation_parameter_t::HUMIDITY, z );
		double T_z = this->getval_( NCPA::attenuation_parameter_t::TEMPERATURE, z );

//		double P_z = this->interpolators_[NCPA::attenuation_parameter_t::PRESSURE]->f(z);
//		double H_z = this->interpolators_[NCPA::attenuation_parameter_t::HUMIDITY]->f(z);
//		double T_z = this->interpolators_[NCPA::attenuation_parameter_t::TEMPERATURE]->f(z);

//		std::cout << "Frequency: " << freq << " Hz" << std::endl
//				<< "Altitude: " << z
//				<< NCPA::Units::toStr(this->get_calculation_units(NCPA::attenuation_parameter_t::ALTITUDE))
//				<< std::endl
//				<< "Pressure: " << P_z
//				<< NCPA::Units::toStr(this->get_calculation_units(NCPA::attenuation_parameter_t::PRESSURE))
//				<< std::endl
//				<< "Temperature: " << T_z
//				<< NCPA::Units::toStr(this->get_calculation_units(NCPA::attenuation_parameter_t::TEMPERATURE))
//				<< std::endl
//				<< "Humidity: " << H_z
//				<< NCPA::Units::toStr(this->get_calculation_units(NCPA::attenuation_parameter_t::HUMIDITY))
//				<< std::endl;
		alpha[i] = NCPA::Units::convert(
				this->NCPA::ISO9613_1AttenuationModel::model( freq, T_z, P_z, H_z ),
				native_output_units_,
				output_units_ );

	}
}

// model assumes temperature in Kelvin, pressure in atmospheres, humidity in percent
double NCPA::ISO9613_1AttenuationModel::model( double freq,
		double temperature, double pressure, double humidity ) {

	double To1 = 273.15;
	double To  = 293.15;

//	double F = freq / pressure;                // frequency per atm

	// calculate saturation pressure
	double Psat = std::pow( 10.0,
		(10.79586*(1-(To1/temperature)) - 5.02808*std::log10(temperature/To1)
			+ 1.50474e-4*(1.0 - std::pow(10.0, -8.29692*((temperature/To1)-1)))
			-  4.2873e-4*(1.0 - std::pow(10.0, -4.76955*((To1/temperature)-1)))
			- 2.2195983) );

	// absolute humidity
	double h = humidity * Psat / pressure;

	// relaxation frequency for nitrogen
	double FrN = pressure * std::sqrt(To/temperature)
		* (9.0 + 280.0*h*std::exp(-4.17*(std::pow(To/temperature,1.0/3.0)-1.0)));

	// relaxation frequency for oxygen
	double FrO = pressure * (24.0 + 4.04e4 * h * (h+0.02)/(h+0.391));

	// attenuation coefficient in nepers/m
	double term1 = 1.84e-11*std::sqrt(temperature/To) / pressure;
	double term2 = 1.275e-2 * std::exp(-2239.1/temperature) / (FrO + freq*freq/FrO);
	double term3 = 1.068e-1 * std::exp(-3352.0 / temperature) / (FrN + freq*freq/FrN);
	return freq * freq * (term1 + std::pow(temperature/To,-2.5)*(term2+term3));
}

NCPA::attenuation_model_t NCPA::ISO9613_1AttenuationModel::type() const {
	return NCPA::attenuation_model_t::ISO9613_1;
}
