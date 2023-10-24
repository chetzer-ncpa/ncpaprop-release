#ifndef NCPA_ATTENUATION_TABULARATTENUATIONMODEL_H_INCLUDED
#define NCPA_ATTENUATION_TABULARATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"
#include "NCPAUnits.h"
#include <vector>


namespace NCPA { class TabularAttenuationModel; }
void swap( NCPA::TabularAttenuationModel &a, NCPA::TabularAttenuationModel &b );

namespace NCPA {
	class TabularAttenuationModel : public AttenuationModel {
	public:
		TabularAttenuationModel() = default;
		TabularAttenuationModel( size_t n, const double *zvec, const double *avec,
				NCPA::units_t z_units );
		TabularAttenuationModel( const TabularAttenuationModel &other );
		TabularAttenuationModel( TabularAttenuationModel &&other );
		friend void ::swap( TabularAttenuationModel &a, TabularAttenuationModel &b );
		TabularAttenuationModel& operator=( TabularAttenuationModel other );
		virtual AttenuationModel* clone() const;
		virtual ~TabularAttenuationModel();

		virtual AttenuationModel* set_parameter(attenuation_parameter_t param, size_t nz,
				const double *p, NCPA::units_t units);
		virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
				const NCPA::VectorWithUnits &p);
		virtual AttenuationModel* set_parameter(attenuation_parameter_t param, double p,
				NCPA::units_t units);
		virtual AttenuationModel* set_parameter(attenuation_parameter_t param,
				const NCPA::ScalarWithUnits &p);

		virtual double attenuation(double f, double z, NCPA::units_t z_units = units_t::NONE);
		virtual void attenuation(size_t n, double f, double *z, double *alpha,
				NCPA::units_t z_units = units_t::NONE);
		virtual units_t get_calculation_units(
				attenuation_parameter_t param ) const;
		virtual const std::vector<attenuation_parameter_t>
			get_required_parameters() const;
		virtual attenuation_model_t type() const;


	protected:
		const std::vector<attenuation_parameter_t> required_parameters_ = {
				attenuation_parameter_t::ATTENUATION
		};
	};
}







#endif
