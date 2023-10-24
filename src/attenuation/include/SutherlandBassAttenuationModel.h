#ifndef NCPA__SUTHERLANDBASSATTENUATIONMODEL_H_INCLUDED
#define NCPA__SUTHERLANDBASSATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"
#include "MoleFractionCalculator.h"
#include "NCPAUnits.h"
#include <vector>

namespace NCPA { class SutherlandBassAttenuationModel; }
void swap( NCPA::SutherlandBassAttenuationModel &a, NCPA::SutherlandBassAttenuationModel &b );

namespace NCPA {

	class SutherlandBassAttenuationModel : public AttenuationModel {
	public:
		SutherlandBassAttenuationModel();
		SutherlandBassAttenuationModel( const SutherlandBassAttenuationModel &other );
		SutherlandBassAttenuationModel( SutherlandBassAttenuationModel &&other );
		virtual ~SutherlandBassAttenuationModel();
		friend void ::swap( SutherlandBassAttenuationModel &a, SutherlandBassAttenuationModel &b );
		SutherlandBassAttenuationModel& operator=( SutherlandBassAttenuationModel other );
		virtual AttenuationModel* clone() const;

		virtual double attenuation(double f, double z, NCPA::units_t z_units = units_t::NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units = units_t::NONE);
		virtual units_t get_calculation_units( attenuation_parameter_t param ) const;
		virtual const std::vector<attenuation_parameter_t> get_required_parameters() const;
		virtual attenuation_model_t type() const;

		static double gamma( double T );
		static double molecular_weight_of_air( double z );
		static double model(
				double freq, double T_z, double P_z, double c_snd_z,
				double X_O2, double X_N2, double X_CO2, double X_O3,
				double X_O, double X_N, double X_H2O );

	protected:

		constexpr static double mu_o  		= 18.192E-6;    // Reference viscosity [kg/(m*s)]
		constexpr static double T_o   		= 293.15;         		// Reference temperature [K]
		constexpr static double P_o   		= 101325;        		// Reference pressure [Pa]
		constexpr static double S     		= 117.0;	     		// Sutherland constant [K]
		constexpr static double R_o		= 8314.48;				// Universal gas constant
//		const double gamma		= 1.4;

		// heat capacity|volume for O2, N2, CO2, and O3
		constexpr static double Cv_R[4] = { 5.0/2.0, 5.0/2.0, 3.0, 3.0 };
		// heat capacity|pressure for O2, N2, CO2, and O3
		constexpr static double Cp_R[4] = { 7.0/2.0, 7.0/2.0, 4.0, 4.0 };
		// characteristic temperature for O2, N2, CO2, and O3
		constexpr static double theta[4] = { 2239.1, 3352.0, 915.0, 1037.0 };

		static constexpr attenuation_parameter_t required_parameters_[2] = {
				attenuation_parameter_t::TEMPERATURE,
				attenuation_parameter_t::PRESSURE
		};

	}; // class SutherlandBassAttenuationModel

} // namespace

#endif
