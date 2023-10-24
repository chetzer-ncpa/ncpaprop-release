/**
 * Derived from:
 *
 * International Organization for Standardization.  1993.  Acoustics -
 *   Attenuation of sound during propagation outdoors - Part 1: Calculation
 *   of the absorption of sound by the atmosphere (ISO Standard No.
 *   9613-1:1993(E)).
 *
 */



#ifndef NCPA_ATTENUATION_ISO96131ATTENUATIONMODEL_H_INCLUDED
#define NCPA_ATTENUATION_ISO96131ATTENUATIONMODEL_H_INCLUDED

#include "NCPAUnits.h"
#include "AttenuationModel.h"
#include <vector>

namespace NCPA { class ISO9613_1AttenuationModel; }
void swap( NCPA::ISO9613_1AttenuationModel &a, NCPA::ISO9613_1AttenuationModel &b );

namespace NCPA {
	class ISO9613_1AttenuationModel : public NCPA::AttenuationModel {
	public:
		ISO9613_1AttenuationModel();
		ISO9613_1AttenuationModel( const ISO9613_1AttenuationModel &other );
		ISO9613_1AttenuationModel( ISO9613_1AttenuationModel &&other );
		virtual ~ISO9613_1AttenuationModel();
		friend void ::swap( ISO9613_1AttenuationModel &a, ISO9613_1AttenuationModel &b );
		ISO9613_1AttenuationModel& operator=( ISO9613_1AttenuationModel other );
		virtual AttenuationModel* clone() const;

		virtual double attenuation(double f, double z, NCPA::units_t z_units = units_t::NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units = units_t::NONE);
		virtual units_t get_calculation_units( attenuation_parameter_t param ) const;
		virtual const std::vector<attenuation_parameter_t> get_required_parameters() const;
		virtual attenuation_model_t type() const;

		// basic model method assumes T in K, P in atm
		static double model( double freq, double temperature, double pressure,
				double humidity );

	protected:
		const std::vector<attenuation_parameter_t> required_parameters_ = {
				attenuation_parameter_t::PRESSURE,
				attenuation_parameter_t::HUMIDITY,
				attenuation_parameter_t::TEMPERATURE
		};
		double getval_( attenuation_parameter_t param, double z );

	};
}









#endif
