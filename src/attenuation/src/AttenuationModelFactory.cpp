#include "AttenuationModelFactory.h"
#include "SutherlandBassAttenuationModel.h"
#include "LosslessAttenuationModel.h"
#include "ISO9613_1AttenuationModel.h"
//#include "CustomAttenuationModel.h"

NCPA::AttenuationModel* NCPA::AttenuationModelFactory::create( attenuation_model_t t ) {
	switch (t) {
	case ATTENUATION_MODEL_NONE:
	case ATTENUATION_MODEL_LOSSLESS:
		return new NCPA::LosslessAttenuationModel();
		break;
	case ATTENUATION_MODEL_SUTHERLAND_BASS:
		return new NCPA::SutherlandBassAttenuationModel();
		break;
	case ATTENUATION_MODEL_CUSTOM:
//		return new NCPA::CustomAttenuationModel();
		return NULL;
		break;
	case ATTENUATION_MODEL_ISO9613_1:
		return new NCPA::ISO9613_1AttenuationModel();
		break;
	default:
		return NULL;
	}
}
