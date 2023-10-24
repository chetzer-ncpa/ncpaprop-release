#ifndef NCPA_ATTENUATION_LOSSLESSATTENUATIONMODEL_H_INCLUDED
#define NCPA_ATTENUATION_LOSSLESSATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"
#include <vector>


namespace NCPA { class LosslessAttenuationModel; }
void swap( NCPA::LosslessAttenuationModel &a, NCPA::LosslessAttenuationModel &b );

namespace NCPA {
	class LosslessAttenuationModel : public AttenuationModel {
	public:
		LosslessAttenuationModel();
		LosslessAttenuationModel( const LosslessAttenuationModel &other );
		LosslessAttenuationModel( LosslessAttenuationModel &&other );
		friend void ::swap( LosslessAttenuationModel &a, LosslessAttenuationModel &b );
		LosslessAttenuationModel& operator=( LosslessAttenuationModel other );
		virtual AttenuationModel* clone() const;
		virtual ~LosslessAttenuationModel();

		virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
				const NCPA::VectorWithUnits &p);

		virtual double attenuation(double f, double z, NCPA::units_t z_units = units_t::NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units = units_t::NONE);
		virtual units_t get_calculation_units(
				attenuation_parameter_t param ) const;
		virtual const std::vector<attenuation_parameter_t>
			get_required_parameters() const;
		virtual attenuation_model_t type() const;

		virtual AttenuationModel* finalize();

	protected:
		const std::vector<attenuation_parameter_t> required_parameters_ = {};
	};
}







#endif
