#include "SutherlandBassAttenuationModel.h"
#include "MoleFractionCalculator.h"
#include "NCPAUnits.h"
#include <cmath>
#include <vector>

NCPA::SutherlandBassAttenuationModel::SutherlandBassAttenuationModel()
		: NCPA::AttenuationModel() {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
}

NCPA::SutherlandBassAttenuationModel::SutherlandBassAttenuationModel(
		const NCPA::SutherlandBassAttenuationModel &other )
		: NCPA::AttenuationModel( other ) {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
}

NCPA::SutherlandBassAttenuationModel::SutherlandBassAttenuationModel(
		NCPA::SutherlandBassAttenuationModel &&other )
		: NCPA::AttenuationModel() {
	native_output_units_ = NCPA::units_t::ATTENUATION_NEPERS_PER_METER;
	::swap( *this, other );
}

void swap( NCPA::SutherlandBassAttenuationModel &a, NCPA::SutherlandBassAttenuationModel &b ) {
	using std::swap;
	::swap( static_cast<NCPA::AttenuationModel&>(a),
			static_cast<NCPA::AttenuationModel&>(b) );
}

NCPA::SutherlandBassAttenuationModel& NCPA::SutherlandBassAttenuationModel::operator=(
		NCPA::SutherlandBassAttenuationModel other ) {
	::swap(*this,other);
	return *this;
}

NCPA::AttenuationModel* NCPA::SutherlandBassAttenuationModel::clone() const {
	return static_cast<NCPA::AttenuationModel *>( new NCPA::SutherlandBassAttenuationModel( *this ) );
}

NCPA::SutherlandBassAttenuationModel::~SutherlandBassAttenuationModel() {}


NCPA::units_t NCPA::SutherlandBassAttenuationModel::get_calculation_units(
		attenuation_parameter_t param ) const {
	switch (param) {
	case NCPA::attenuation_parameter_t::ALTITUDE:		// default km
		return NCPA::units_t::DISTANCE_KILOMETERS;
		break;
	case NCPA::attenuation_parameter_t::TEMPERATURE:	// default K
		return NCPA::units_t::TEMPERATURE_KELVIN;
		break;
	case NCPA::attenuation_parameter_t::PRESSURE:		// default Pa
		return NCPA::units_t::PRESSURE_PASCALS;
		break;
	default:
		return NCPA::units_t::NONE;
	}
}

const std::vector<NCPA::attenuation_parameter_t>
NCPA::SutherlandBassAttenuationModel::get_required_parameters() const {
	return std::vector<NCPA::attenuation_parameter_t>(
			NCPA::SutherlandBassAttenuationModel::required_parameters_,
			NCPA::SutherlandBassAttenuationModel::required_parameters_ + 2);
}

double NCPA::SutherlandBassAttenuationModel::attenuation(
		double f, double z, NCPA::units_t z_units) {
	double alpha;
	this->attenuation( 1, f, &z, &alpha, z_units );
	return alpha;
}

void NCPA::SutherlandBassAttenuationModel::attenuation(
		size_t n, double freq, double *zvec, double *alpha, NCPA::units_t z_units) {
	if (z_units == NCPA::units_t::NONE) {
		z_units = this->get_calculation_units(NCPA::attenuation_parameter_t::ALTITUDE);
	}

	std::map<NCPA::molecular_constituent_t,NCPA::MoleFractionCalculator*> frac_map
		= NCPA::MoleFractionCalculator::build_all();

	for (size_t i = 0; i < n; i++) {
		double z = NCPA::Units::convert(zvec[i],z_units,
				this->get_calculation_units(NCPA::attenuation_parameter_t::ALTITUDE) );
		double P_z = this->interpolators_[NCPA::attenuation_parameter_t::PRESSURE]->f(z);
		double T_z = this->interpolators_[NCPA::attenuation_parameter_t::TEMPERATURE]->f(z);

		// calculated parameters
		double c_snd_z;
		auto density_it = this->interpolators_.find(NCPA::attenuation_parameter_t::DENSITY);
		if (density_it != this->interpolators_.end()) {
			c_snd_z = std::sqrt( NCPA::SutherlandBassAttenuationModel::gamma(T_z) * P_z
					/ density_it->second->f(z) );  // m/s
		} else {
			c_snd_z = std::sqrt(NCPA::SutherlandBassAttenuationModel::gamma(T_z) * R_o * T_z
					/ NCPA::SutherlandBassAttenuationModel::molecular_weight_of_air(z) );
		}

		// gas fractions
		std::map<NCPA::molecular_constituent_t,double> gasfractions;
		for (auto fit = frac_map.begin(); fit != frac_map.end(); ++fit) {
			gasfractions[fit->first] = fit->second->calculate(z,
					this->get_calculation_units(NCPA::attenuation_parameter_t::ALTITUDE) );
		}

//		double mu      = mu_o * std::sqrt(
//			T_z/T_o ) * ( (1.0 + S/T_o) / (1.0 + S/T_z)
//			); // Viscosity [kg/(m*s)]
//		double nu      = ( 8.0 * M_PI * freq * mu ) / ( 3.0 * P_z ); // Nondimensional frequency

		// Gas fractions
//		double X[7];
//
//		//-------- Gas fraction polynomial fits -----------------------------------
//		if (z > 90.0) {                                         // O2 profile
//			X[0] = std::pow( 10.0,
//						49.296 - (1.5524*z) + (1.8714E-2*std::pow(z,2))
//				   		- (1.1069E-4*std::pow(z,3)) + (3.199E-7*std::pow(z,4))
//				   		- (3.6211E-10*std::pow(z,5))
//				);
//		} else {
//			X[0] = std::pow(10.0,-0.67887);
//		}
//
//		if (z > 76.0) {                                        // N2 profile
//			X[1] = std::pow(10.0,(1.3972E-1)-(5.6269E-3*z) + (3.9407E-5*std::pow(z,2) )
//				   - (1.0737E-7*std::pow(z,3)));
//		} else {
//			X[1] = std::pow(10.0,-0.10744);
//		}
//
//		X[2] = std::pow(10.0,-3.3979);                              // CO2 profile
//
//		if (z > 80.0) {                                        // O3 profile
//			X[3] = std::pow(10.0,-4.234-(3.0975E-2*z));
//		} else {
//			X[3] = std::pow(10.0,-19.027+(1.3093*z) - (4.6496E-2*std::pow(z,2))
//				   + (7.8543E-4*std::pow(z,3)) - (6.5169E-6*std::pow(z,4))
//				   + (2.1343E-8*std::pow(z,5)));
//		}
//
//		if (z > 95.0 ) {                                       // O profile
//			X[4] = std::pow(10.0,-3.2456+(4.6642E-2*z)-(2.6894E-4*std::pow(z,2))+(5.264E-7*std::pow(z,3)));
//		} else {
//			X[4] = std::pow(10.0,-11.195+(1.5408E-1*z)-(1.4348E-3*std::pow(z,2))+(1.0166E-5*std::pow(z,3)));
//		}
//
//		// N profile
//		X[5]  = std::pow(10.0,-53.746+(1.5439*z)-(1.8824E-2*std::pow(z,2))+(1.1587E-4*std::pow(z,3))
//			    -(3.5399E-7*std::pow(z,4))+(4.2609E-10*std::pow(z,5)));
//
//		if (z > 30.0 ) {                                        // H2O profile
//			X[6] = std::pow(10.0,-4.2563+(7.6245E-2*z)-(2.1824E-3*std::pow(z,2))-(2.3010E-6*std::pow(z,3))
//				   +(2.4265E-7*std::pow(z,4))-(1.2500E-09*std::pow(z,5)));
//		} else {
//			if (z > 100.0) {
//				X[6] = std::pow(10.0,-0.62534-(8.3665E-2*z));
//			} else {
//				X[6] = std::pow(10.0,-1.7491+(4.4986E-2*z)-(6.8549E-2*std::pow(z,2))
//					   +(5.4639E-3*std::pow(z,3))-(1.5539E-4*std::pow(z,4))
//					   +(1.5063E-06*std::pow(z,5)));
//			}
//		}

		alpha[i] = NCPA::SutherlandBassAttenuationModel::model(
				z, T_z, P_z, c_snd_z,
				gasfractions[NCPA::molecular_constituent_t::O2],
				gasfractions[NCPA::molecular_constituent_t::N2],
				gasfractions[NCPA::molecular_constituent_t::CO2],
				gasfractions[NCPA::molecular_constituent_t::O3],
				gasfractions[NCPA::molecular_constituent_t::O],
				gasfractions[NCPA::molecular_constituent_t::N],
				gasfractions[NCPA::molecular_constituent_t::H2O]
		);
	}
}

double NCPA::SutherlandBassAttenuationModel::gamma( double T ) {
	const double A[6] = { 1.371, 2.460e-4, -6.436e-7, 5.2e-10, -1.796e-13, 2.182e-17 };
	return NCPA::evalpoly( 6, A, T );
}

double NCPA::SutherlandBassAttenuationModel::molecular_weight_of_air( double z ) {
	std::vector<NCPA::molecular_constituent_t> v = NCPA::MoleFractionCalculator::available_types();
	double num = 0.0, denom = 0.0;
	for (auto it = v.begin(); it != v.end(); ++it) {
		NCPA::MoleFractionCalculator *calc = NCPA::MoleFractionCalculator::build(*it);
		double Mi = calc->molecular_weight();
		double Xi = calc->calculate(z, NCPA::units_t::DISTANCE_KILOMETERS);
		num += Xi*Mi;
		denom += Xi;
	}
	return num/denom;
}

NCPA::attenuation_model_t NCPA::SutherlandBassAttenuationModel::type() const {
	return NCPA::attenuation_model_t::SUTHERLAND_BASS;
}

double NCPA::SutherlandBassAttenuationModel::model(
		double freq, double T_z, double P_z, double c_snd_z,
		double X_O2, double X_N2, double X_CO2, double X_O3,
		double X_O, double X_N, double X_H2O ) {
	double mu      = mu_o * std::sqrt(
		T_z/T_o ) * ( (1.0 + S/T_o) / (1.0 + S/T_z)
	); // Viscosity [kg/(m*s)], Eq. 17

	double nu      = ( 8.0 * M_PI * freq * mu ) / ( 3.0 * P_z ); // Nondimensional frequency

	double X_ON = (X_O2 + X_N2)/0.9903;

	//-------- Rotational collision number-------------------------------------
	double Z_rot_0 = 54.1*std::exp(-17.3*(std::pow(T_z,-1.0/3.0)));   // O2
	double Z_rot_1 = 63.3*std::exp(-16.7*(std::pow(T_z,-1.0/3.0)));   // N2
	double Z_rot_   = 1.0 / ( (X_N2 / Z_rot_1) + (X_O2 / Z_rot_0) );	// Eq. 13

	//-------- Nondimensional atmospheric quantities---------------------------
	double sigma = 5.0 / std::sqrt(21.0);
	double nn    = (4.0/5.0) * std::sqrt(3.0/7.0) * Z_rot_;
	double chi   = 3.0 * nn * nu / 4.0;
	double cchi  = 2.36 * chi;

	//---------Classical + rotational loss/dispersion--------------------------
	double a_cl =
			(2 * M_PI * freq / c_snd_z)
			* std::sqrt(
					0.5
					* (std::sqrt(1+std::pow(nu,2))-1)
					* (1+std::pow(cchi,2))
					/ (
							(1+std::pow(nu,2))
							* (1+std::pow(sigma*cchi,2))
					)
			);  // Eq. 16
	double a_rot   =
			(2*M_PI*freq/c_snd_z)
			* X_ON
			* (
					(pow(sigma,2)-1)
					* chi
					/ (2*sigma)
			)
			* sqrt(
					0.5
					* (sqrt(1+pow(nu,2))+1)
					/(
							(1+pow(nu,2))
							* (1+pow(sigma*cchi,2))
							* (1+pow(cchi,2))
					)
			);  // Eq. 18
	double a_diff  = 0.003*a_cl;	// Eq. 25

	//---------Vibrational relaxation-------------------------------------------
	double Tr = std::pow( T_z/T_o ,-1.0/3.0 ) - 1.0;
	double A1 = (X_O2+X_N2) * 24.0 * std::exp(-9.16*Tr);
	double A2 = (X_O+X_N) * 2400.0;
	double B  = 40400.0 * std::exp( 10.0*Tr );
	double C  = 0.02 * std::exp( -11.2*Tr );
	double D  = 0.391 * std::exp( 8.41*Tr );
	double E  = 9.0 * std::exp( -19.9*Tr );
	double F  = 60000.0;
	double G  = 28000.0 * std::exp( -4.17*Tr );
	double H  = 22000.0 * std::exp( -7.68*Tr );
	double I  = 15100.0 * std::exp( -10.4*Tr );
	double J  = 11500.0 * std::exp( -9.17*Tr );
	double K  = (8.48E08) * std::exp( 9.17*Tr );
	double L  = std::exp( -7.72*Tr );
	double ZZ = H*X_CO2 + I*(X_O2 + 0.5*X_O) + J*(X_N2 + 0.5*X_N ) + K*(X_H2O + X_O3);
	double hu = 100.0 * (X_O3 + X_H2O);
	double f_vib[ 4 ], a_vib_c[ 4 ];
	f_vib[0] = (P_z/P_o) * (mu_o/mu) * (A1 + A2 + B*hu*(C+hu)*(D+hu));	// Eq. 28
	f_vib[1] = (P_z/P_o) * (mu_o/mu) * (E + F*X_O3 + G*X_H2O);			// Eq. 29
	f_vib[2] = (P_z/P_o) * (mu_o/mu) * ZZ;								// Eq. 30
	f_vib[3] = (P_z/P_o) * (mu_o/mu) * (1.2E5)*L;

	double a_vib = 0.0;
	double X[4] = { X_O2, X_N2, X_CO2, X_O3 };
	for (size_t m=0; m<4; m++)
	{
		double thetaT = theta[m] / T_z;
		double C_R =
				std::pow(thetaT,2.0)
				* std::exp(-thetaT)
				/ std::pow(1 - std::exp(-thetaT), 2.0);	// Eq. 27
		double A_max =
				X[m] * (M_PI/2.0) * C_R
				/ ( Cp_R[m] * (Cv_R[m] + C_R) );	// Eq. 26
		a_vib_c[m] =
				(A_max/c_snd_z)
				* (
						(2*(pow(freq,2.0)) / f_vib[m])
						/ (1 + pow(freq/f_vib[m],2.0)));
		a_vib += a_vib_c[m];	// Eq. 1
	}

	return a_cl + a_rot + a_diff + a_vib;
}
