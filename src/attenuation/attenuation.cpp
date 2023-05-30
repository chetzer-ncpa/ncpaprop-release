#include "attenuation.h"
#include <vector>
#include "units.h"
#include "VectorWithUnits.h"

// Base class methods
NCPA::AttenuationModel::~AttenuationModel() {

}

void NCPA::AttenuationModel::set_vertical_scale(size_t nz, const double *z_km) {
	NCPA::VectorWithUnits *zvec = new NCPA::VectorWithUnits(nz, z_km, NCPA_ATTENUATION_DEFAULT_UNITS_Z );
	parameter_map.insert( attenuation_parameter_pair_t(ATTENUATION_MODEL_PARAMETER_ALTITUDE,zvec) );
}
