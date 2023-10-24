#include <vector>
#include <stdexcept>

#include "AttenuationModel.h"
#include "LosslessAttenuationModel.h"
#include "ISO9613_1AttenuationModel.h"
#include "SutherlandBassAttenuationModel.h"
#include "TabularAttenuationModel.h"
#include "NCPAUnits.h"

std::string NCPA::AttenuationModel::as_string(NCPA::attenuation_parameter_t t) {
	switch (t) {
	case NCPA::attenuation_parameter_t::ALTITUDE:
		return "altitude";
		break;
	case NCPA::attenuation_parameter_t::TEMPERATURE:
		return "temperature";
		break;
	case NCPA::attenuation_parameter_t::DENSITY:
		return "density";
		break;
	case NCPA::attenuation_parameter_t::HUMIDITY:
		return "humidity";
		break;
	case NCPA::attenuation_parameter_t::PRESSURE:
		return "pressure";
		break;
	case NCPA::attenuation_parameter_t::ATTENUATION:
		return "attenuation";
		break;
	default:
		throw std::out_of_range( "Unrecognized parameter type" );
	}
}

// Base class methods
NCPA::AttenuationModel::AttenuationModel( const NCPA::AttenuationModel &other ) {
	this->vector_parameters_ = other.vector_parameters_;
	this->scalar_parameters_ = other.scalar_parameters_;
	for (auto it = other.interpolators_.cbegin(); it != other.interpolators_.cend(); ++it) {
		this->interpolators_[ it->first ] = it->second->clone();
	}
	this->z_vector_ = other.z_vector_;
	this->interpolator_type_ = other.interpolator_type_;
	this->finalized_ = other.finalized_;
}

NCPA::AttenuationModel::AttenuationModel( NCPA::AttenuationModel &&other ) {
	::swap(*this,other);
}

void swap( NCPA::AttenuationModel &a, NCPA::AttenuationModel &b ) {
	using std::swap;
	swap( a.vector_parameters_, b.vector_parameters_ );
	swap( a.scalar_parameters_, b.scalar_parameters_ );
	swap( a.interpolators_, b.interpolators_ );
	swap( a.z_vector_, b.z_vector_ );
	swap( a.interpolator_type_, b.interpolator_type_ );
	swap( a.finalized_, b.finalized_ );
}

NCPA::AttenuationModel::~AttenuationModel() {}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_altitude_parameter(const NCPA::VectorWithUnits &z) {
	this->set_parameter_(NCPA::attenuation_parameter_t::ALTITUDE, z);
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_altitude_parameter(size_t nz, const double *p,
		NCPA::units_t units) {
	return this->set_altitude_parameter( NCPA::VectorWithUnits( nz, p, units ) );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_altitude_parameter(const std::vector<double> &p,
		NCPA::units_t units) {
	return this->set_altitude_parameter( NCPA::VectorWithUnits( p, units ) );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_altitude_parameter(size_t nz, const double *p,
		const std::string &units) {
	return this->set_altitude_parameter( NCPA::VectorWithUnits( nz, p, units ) );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_altitude_parameter(const std::vector<double> &p,
		const std::string &units) {
	return this->set_altitude_parameter( NCPA::VectorWithUnits( p, units ) );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		const NCPA::VectorWithUnits &p) {
	if (z_vector_.empty()) {
		throw std::out_of_range( "AttenuationModel: No existing vertical scale set, so can't interpret new vertical parameter." );
	}
	if (z_vector_.size() != p.size()) {
		throw std::out_of_range( "AttenuationModel: Existing vertical scale has different number of points than new vertical parameter." );
	}
	this->set_parameter_(param, p);
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		size_t nz, const double *p, NCPA::units_t units) {
	return this->set_parameter(param, NCPA::VectorWithUnits( nz, p, units ));
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		size_t nz, const double *p, const std::string &units) {
	return this->set_parameter(param, NCPA::VectorWithUnits( nz, p, units ));
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		const std::vector<double> &p, NCPA::units_t units) {
	return this->set_parameter(param, NCPA::VectorWithUnits( p, units ));
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		const std::vector<double> &p, const std::string &units) {
	return this->set_parameter(param, NCPA::VectorWithUnits( p, units ));
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		const NCPA::ScalarWithUnits &p) {
	this->set_parameter_(param, p);
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		double p, NCPA::units_t units) {
	return this->set_parameter(param, NCPA::ScalarWithUnits( p, units ));
}

NCPA::AttenuationModel* NCPA::AttenuationModel::set_parameter(attenuation_parameter_t param,
		double p, const std::string &units) {
	return this->set_parameter(param, NCPA::ScalarWithUnits( p, units ));
}


void NCPA::AttenuationModel::set_parameter_(
		attenuation_parameter_t param,
		NCPA::VectorWithUnits p) {
	if (finalized_) {
		throw std::logic_error("Attenuation model has already been finalized, cannot modify!");
	}
	if (param == NCPA::attenuation_parameter_t::ALTITUDE) {
		z_vector_ = p;
	} else {
		remove_parameter_( param );
		vector_parameters_[ param ] = p;
	}
}

void NCPA::AttenuationModel::set_parameter_(
		attenuation_parameter_t param,
		NCPA::ScalarWithUnits p) {
	if (finalized_) {
		throw std::logic_error("Attenuation model has already been finalized, cannot modify!");
	}
	remove_parameter_( param );
	scalar_parameters_[ param ] = p;
}

void NCPA::AttenuationModel::remove_parameter_( attenuation_parameter_t param ) {
	if (finalized_) {
		throw std::logic_error("Attenuation model has already been finalized, cannot modify!");
	}
	if (vector_parameters_.find(param) != vector_parameters_.end()) {
		vector_parameters_.erase(param);
	}
	if (scalar_parameters_.find(param) != scalar_parameters_.end()) {
		scalar_parameters_.erase(param);
	}
}

NCPA::AttenuationModel* NCPA::AttenuationModel::finalize() {
	if (z_vector_.size() == 0) {
		throw std::logic_error( "No altitude vector has been set!" );
	}

	std::vector<NCPA::attenuation_parameter_t> reqs = this->get_required_parameters();
	for (auto it = reqs.cbegin(); it != reqs.cend(); ++it) {
		if (this->vector_parameters_.find(*it) == this->vector_parameters_.end()
				&& this->scalar_parameters_.find(*it) == this->scalar_parameters_.end()) {
			throw std::logic_error( "Missing parameter: " + NCPA::AttenuationModel::as_string(*it) );
		}
	}

	clear_interpolators_();
	z_vector_.convert_units(
			this->get_calculation_units( NCPA::attenuation_parameter_t::ALTITUDE ) );
	for (auto it = this->vector_parameters_.begin(); it != this->vector_parameters_.end(); ++it) {
		it->second.convert_units( this->get_calculation_units( it->first ) );
		this->interpolators_[ it->first ] = NCPA::Interpolator1D::build(this->interpolator());
		this->interpolators_[ it->first ]->init();
		this->interpolators_[ it->first ]->allocate(z_vector_.size());
		this->interpolators_[ it->first ]->set( z_vector_.as_doubles(), it->second.as_doubles() );
		this->interpolators_[ it->first ]->ready();
	}
	for (auto it = this->scalar_parameters_.begin(); it != this->scalar_parameters_.end(); ++it) {
		it->second.convert_units( this->get_calculation_units( it->first ) );
	}
	this->finalized_ = true;
	return static_cast<NCPA::AttenuationModel *>( this );
}

void NCPA::AttenuationModel::clear_interpolators_() {
	for (auto it = interpolators_.begin(); it != interpolators_.end(); ++it) {
		if (it->second != nullptr) {
			delete it->second;
			it->second = nullptr;
		}
	}
	interpolators_.clear();
}

std::string NCPA::AttenuationModel::as_string( NCPA::attenuation_model_t t ) {
	switch (t) {
	case NCPA::attenuation_model_t::LOSSLESS:
		return "Lossless attenuation model";
		break;
	case NCPA::attenuation_model_t::ISO9613_1:
		return "ISO 9613-1 attenuation model";
		break;
	case NCPA::attenuation_model_t::SUTHERLAND_BASS:
		return "Sutherland/Bass attenuation model";
		break;
	case NCPA::attenuation_model_t::TABULAR:
		return "Tabular attenuation model";
		break;
	default:
		throw std::out_of_range("Unrecognized, invalid or unimplemented.");
	}
}

bool NCPA::AttenuationModel::can_build( NCPA::attenuation_model_t t ) {
	switch (t) {
	case NCPA::attenuation_model_t::LOSSLESS:
	case NCPA::attenuation_model_t::ISO9613_1:
	case NCPA::attenuation_model_t::SUTHERLAND_BASS:
	case NCPA::attenuation_model_t::TABULAR:
		return true;
		break;
	default:
		return false;
	}
}

NCPA::AttenuationModel *NCPA::AttenuationModel::build( NCPA::attenuation_model_t t ) {
	switch (t) {
	case NCPA::attenuation_model_t::LOSSLESS:
		return static_cast<NCPA::AttenuationModel*>(new NCPA::LosslessAttenuationModel());
		break;
	case NCPA::attenuation_model_t::ISO9613_1:
		return static_cast<NCPA::AttenuationModel*>(new NCPA::ISO9613_1AttenuationModel());
		break;
	case NCPA::attenuation_model_t::SUTHERLAND_BASS:
		return static_cast<NCPA::AttenuationModel*>(new NCPA::SutherlandBassAttenuationModel());
		break;
	case NCPA::attenuation_model_t::TABULAR:
		return static_cast<NCPA::AttenuationModel*>(new NCPA::TabularAttenuationModel());
		break;
	default:
		throw std::out_of_range("Unrecognized or unimplemented interpolator requested.");
	}
}

NCPA::units_t NCPA::AttenuationModel::output_units() const {
	return output_units_ == NCPA::units_t::NONE ? native_output_units_ : output_units_;
}

NCPA::AttenuationModel* NCPA::AttenuationModel::output_units( NCPA::units_t u ) {
	output_units_ = u;
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::AttenuationModel* NCPA::AttenuationModel::interpolator( NCPA::interpolator1d_t inttype ) {
	interpolator_type_ = inttype;
	return static_cast<NCPA::AttenuationModel *>( this );
}

NCPA::interpolator1d_t NCPA::AttenuationModel::interpolator() const {
	return interpolator_type_;
}





