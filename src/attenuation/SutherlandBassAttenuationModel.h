#ifndef NCPA__SUTHERLANDBASSATTENUATIONMODEL_H_INCLUDED
#define NCPA__SUTHERLANDBASSATTENUATIONMODEL_H_INCLUDED

#include "AttenuationModel.h"

namespace NCPA {

	class SutherlandBassAttenuationModel : public AttenuationModel {
	public:
		SutherlandBassAttenuationModel();
		virtual ~SutherlandBassAttenuationModel();

		virtual double attenuation(double z, double f);
		virtual double attenuation(double z, double f,
				std::map<attenuation_model_parameter_t,double> params);
	}; // class SutherlandBassAttenuationModel

} // namespace

#endif
