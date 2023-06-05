#include "SutherlandBassAttenuationModel.h"
#include "NCPAUnits.h"
#include <cmath>
#include <vector>


NCPA::units_t NCPA::SutherlandBassAttenuationModel::get_calculation_units(
		attenuation_model_parameter_label_t param ) const {
	switch (param) {
	case ATTENUATION_MODEL_PARAMETER_ALTITUDE:		// default km
		return NCPA::UNITS_DISTANCE_KILOMETERS;
		break;
	case ATTENUATION_MODEL_PARAMETER_TEMPERATURE:	// default K
		return NCPA::UNITS_TEMPERATURE_KELVIN;
		break;
	case ATTENUATION_MODEL_PARAMETER_DENSITY:		// default kg/m3
		return NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER;
		break;
	case ATTENUATION_MODEL_PARAMETER_PRESSURE:		// default Pa
		return NCPA::UNITS_PRESSURE_PASCALS;
		break;
	default:
		return NCPA::UNITS_NONE;
	}
}

const std::vector<NCPA::attenuation_model_parameter_label_t>
NCPA::SutherlandBassAttenuationModel::get_required_parameters() const {
	return required_parameters;
}

double NCPA::SutherlandBassAttenuationModel::attenuation(
		double f, double z, NCPA::units_t z_units) {
	double alpha;
	this->attenuation( 1, f, &z, &alpha, z_units );
	return alpha;
}

void NCPA::SutherlandBassAttenuationModel::attenuation(
		size_t n, double f, double *zvec, double *alpha, NCPA::units_t z_units) {
	if (z_units == NCPA::UNITS_NONE) {
		z_units = this->get_calculation_units(ATTENUATION_MODEL_PARAMETER_ALTITUDE);
	}

	// get units right
	for (auto it = get_required_parameters().cbegin(); it != get_required_parameters().cend(); ++it) {
		this->parameter_map[*it]->convert_units( this->get_calculation_units( *it ) );
	}

	for (size_t i = 0; i < n; i++) {
		double z = NCPA::Units::convert(zvec[i],z_units,
				this->get_calculation_units(ATTENUATION_MODEL_PARAMETER_ALTITUDE) );
		double P_z = this->parameter_map[ATTENUATION_MODEL_PARAMETER_PRESSURE]->get(z);
		double D_z = this->parameter_map[ATTENUATION_MODEL_PARAMETER_DENSITY]->get(z);
		double T_z = this->parameter_map[ATTENUATION_MODEL_PARAMETER_TEMPERATURE]->get(z);

		// calculated parameters
		double c_snd_z = std::sqrt( gamma * P_z / D_z );  // m/s
		double mu      = mu_o * std::sqrt(
			T_z/T_o ) * ( (1.0 + S/T_o) / (1.0 + S/T_z)
			); // Viscosity [kg/(m*s)]
		double nu      = ( 8.0 * PI * freq * mu ) / ( 3.0 * P_z ); // Nondimensional frequency

		// Gas fractions
		double X[7];

		//-------- Gas fraction polynomial fits -----------------------------------
		if (z > 90.0) {                                         // O2 profile
			X[0] = std::pow( 10.0,
						49.296 - (1.5524*z) + (1.8714E-2*std::pow(z,2))
				   		- (1.1069E-4*std::pow(z,3)) + (3.199E-7*std::pow(z,4))
				   		- (3.6211E-10*std::pow(z,5))
				);
		} else {
			X[0] = std::pow(10.0,-0.67887);
		}

		if (z > 76.0) {                                        // N2 profile
			X[1] = std::pow(10.0,(1.3972E-1)-(5.6269E-3*z) + (3.9407E-5*std::pow(z,2) )
				   - (1.0737E-7*std::pow(z,3)));
		} else {
			X[1] = std::pow(10.0,-0.10744);
		}

		X[2] = std::pow(10.0,-3.3979);                              // CO2 profile

		if (z > 80.0) {                                        // O3 profile
			X[3] = std::pow(10.0,-4.234-(3.0975E-2*z));
		} else {
			X[3] = std::pow(10.0,-19.027+(1.3093*z) - (4.6496E-2*std::pow(z,2))
				   + (7.8543E-4*std::pow(z,3)) - (6.5169E-6*std::pow(z,4))
				   + (2.1343E-8*std::pow(z,5)));
		}

		if (z > 95.0 ) {                                       // O profile
			X[4] = std::pow(10.0,-3.2456+(4.6642E-2*z)-(2.6894E-4*std::pow(z,2))+(5.264E-7*std::pow(z,3)));
		} else {
			X[4] = std::pow(10.0,-11.195+(1.5408E-1*z)-(1.4348E-3*std::pow(z,2))+(1.0166E-5*std::pow(z,3)));
		}

		// N profile
		X[5]  = std::pow(10.0,-53.746+(1.5439*z)-(1.8824E-2*std::pow(z,2))+(1.1587E-4*std::pow(z,3))
			    -(3.5399E-7*std::pow(z,4))+(4.2609E-10*std::pow(z,5)));

		if (z > 30.0 ) {                                        // H2O profile
			X[6] = std::pow(10.0,-4.2563+(7.6245E-2*z)-(2.1824E-3*std::pow(z,2))-(2.3010E-6*std::pow(z,3))
				   +(2.4265E-7*std::pow(z,4))-(1.2500E-09*std::pow(z,5)));
		} else {
			if (z > 100.0) {
				X[6] = std::pow(10.0,-0.62534-(8.3665E-2*z));
			} else {
				X[6] = std::pow(10.0,-1.7491+(4.4986E-2*z)-(6.8549E-2*std::pow(z,2))
					   +(5.4639E-3*std::pow(z,3))-(1.5539E-4*std::pow(z,4))
					   +(1.5063E-06*std::pow(z,5)));
			}
		}
		double X_ON = (X[0] + X[1])/0.9903;

		//-------- Rotational collision number-------------------------------------
		double Z_rot_0 = 54.1*std::exp(-17.3*(std::pow(T_z,-1.0/3.0)));   // O2
		double Z_rot_1 = 63.3*std::exp(-16.7*(std::pow(T_z,-1.0/3.0)));   // N2
		double Z_rot_   = 1.0 / ( (X[1] / Z_rot_1) + (X[0] / Z_rot_0) );

		//-------- Nondimensional atmospheric quantities---------------------------
		double sigma = 5.0 / std::sqrt(21.0);
		double nn    = (4.0/5.0) * std::sqrt(3.0/7.0) * Z_rot_;
		double chi   = 3.0 * nn * nu / 4.0;
		double cchi  = 2.36 * chi;

		//---------Classical + rotational loss/dispersion--------------------------
		double a_cl    = (2 * PI * freq / c_snd_z)
							* std::sqrt( 0.5 * (std::sqrt(1+std::pow(nu,2))-1) * (1+std::pow(cchi,2))
							/ ((1+std::pow(nu,2))*(1+std::pow(sigma*cchi,2))) );
		double a_rot   = (2*PI*freq/c_snd_z) * X_ON
							* ((pow(sigma,2)-1)*chi/(2*sigma))
							* sqrt(0.5*(sqrt(1+pow(nu,2))+1)/((1+pow(nu,2))*(1+pow(cchi,2))));
		double a_diff  = 0.003*a_cl;

		//---------Vibrational relaxation-------------------------------------------
		double Tr = std::pow( T_z/T_o ,-1.0/3.0 ) - 1.0;
		double A1 = (X[0]+X[1]) * 24.0 * std::exp(-9.16*Tr);
		double A2 = (X[4]+X[5]) * 2400.0;
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
		double ZZ = H*X[2] + I*(X[0] + 0.5*X[4]) + J*(X[1] + 0.5*X[5] ) + K*(X[6] + X[3]);
		double hu = 100.0 * (X[3] + X[6]);
		double f_vib[ 4 ], a_vib_c[ 4 ];
		f_vib[0] = (P_z/P_o) * (mu_o/mu) * (A1 + A2 + B*hu*(C+hu)*(D+hu));
		f_vib[1] = (P_z/P_o) * (mu_o/mu) * (E + F*X[3] + G*X[6]);
		f_vib[2] = (P_z/P_o) * (mu_o/mu) * ZZ;
		f_vib[3] = (P_z/P_o) * (mu_o/mu) * (1.2E5)*L;

		double a_vib = 0.0;
		for (size_t m=0; m<4; m++)
		{
			double C_R        = ((std::pow(theta[m]/T_z,2))*std::exp(-theta[m]/T_z))/(std::pow(1-std::exp(-theta[m]/T_z),2));
			double A_max      = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
			       // A_max      = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
			a_vib_c[m] = (A_max/c_snd_z)*((2*(pow(freq,2))/f_vib[m])/(1+pow(freq/f_vib[m],2)));
			a_vib      = a_vib + a_vib_c[m];
		}

		alpha[i] = a_cl + a_rot + a_diff + a_vib;
	}

}
