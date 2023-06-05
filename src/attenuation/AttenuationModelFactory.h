#ifndef NCPA__ATTENUATIONMODELFACTORY_H_INCLUDED
#define NCPA__ATTENUATIONMODELFACTORY_H_INCLUDED

#include "AttenuationModel.h"
#include "SutherlandBassAttenuationModel.h"
#include <unordered_map>

namespace NCPA {
	enum attenuation_model_t : unsigned int {
		ATTENUATION_MODEL_NONE = 0,
		ATTENUATION_MODEL_LOSSLESS,	// same as none
		ATTENUATION_MODEL_SUTHERLAND_BASS,
		ATTENUATION_MODEL_CUSTOM,
		ATTENUATION_MODEL_ISO9613_1
	};

	class AttenuationModelFactory {
	public:
		static AttenuationModel *create( attenuation_model_t t );
	};
}

#endif
