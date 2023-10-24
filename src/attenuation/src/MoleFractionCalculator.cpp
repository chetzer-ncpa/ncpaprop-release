#include "MoleFractionCalculator.h"

#include "NCPAUnits.h"
#include "NCPACommon.h"

#include <stdexcept>
#include <cmath>
#include <vector>

NCPA::MoleFractionCalculator::MoleFractionCalculator( double molweight )
	: molecular_weight_(molweight) {}

double NCPA::MoleFractionCalculator::calculate( double z, NCPA::units_t units ) {
	double z_km = NCPA::Units::convert( z, units, NCPA::units_t::DISTANCE_KILOMETERS );
	if (poly_ != NULL) {
		return std::pow( 10.0, poly_->calculate( z_km ) );
	} else {
		throw std::invalid_argument( "Calculator has not been initialized" );
	}
}

NCPA::MoleFractionCalculator::~MoleFractionCalculator() {
	if (poly_ != NULL) {
		delete poly_;
	}
}

double NCPA::MoleFractionCalculator::molecular_weight() const { return molecular_weight_; }

std::vector<NCPA::molecular_constituent_t> NCPA::MoleFractionCalculator::available_types() {
	NCPA::molecular_constituent_t types[8] = {
			NCPA::molecular_constituent_t::O2,
			NCPA::molecular_constituent_t::N2,
			NCPA::molecular_constituent_t::CO2,
			NCPA::molecular_constituent_t::O3,
			NCPA::molecular_constituent_t::O,
			NCPA::molecular_constituent_t::N,
			NCPA::molecular_constituent_t::H2O,
			NCPA::molecular_constituent_t::AR
	};
	return std::vector<NCPA::molecular_constituent_t>( types, types+8 );
}

std::map<NCPA::molecular_constituent_t,NCPA::MoleFractionCalculator*>
NCPA::MoleFractionCalculator::build_all() {
	std::map<NCPA::molecular_constituent_t,NCPA::MoleFractionCalculator*> fracs;
	auto vtypes = NCPA::MoleFractionCalculator::available_types();
	for (auto vit = vtypes.cbegin(); vit != vtypes.cend(); ++vit) {
		fracs[*vit] = NCPA::MoleFractionCalculator::build( *vit );
	}
	return fracs;
}

NCPA::MoleFractionCalculator * NCPA::MoleFractionCalculator::build(NCPA::molecular_constituent_t type) {
	switch (type) {
	case NCPA::molecular_constituent_t::O2:
		return new NCPA::O2MoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::N2:
		return new NCPA::N2MoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::O3:
		return new NCPA::O3MoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::CO2:
		return new NCPA::CO2MoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::O:
		return new NCPA::OMoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::N:
		return new NCPA::NMoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::H2O:
		return new NCPA::H2OMoleFractionCalculator();
		break;
	case NCPA::molecular_constituent_t::AR:
		return new NCPA::ArMoleFractionCalculator();
		break;
	default:
		return NULL;
	}
}


NCPA::O2MoleFractionCalculator::O2MoleFractionCalculator() : NCPA::MoleFractionCalculator(31.9988) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -0.67887 };
	double c2[ 6 ] = { 49.296, -1.5524, 1.8714e-2, -1.1069e-4, 3.199e-7, -3.6211e-10 };
	poly_->assign( -9999.9, 90.0, 1, c1 );
	poly_->assign( 90.0, 240.0, 6, c2 );
}
//static double NCPA::O2MoleFractionCalculator::molecular_weight() const { return 31.9988; }

NCPA::N2MoleFractionCalculator::N2MoleFractionCalculator() : NCPA::MoleFractionCalculator(28.0134) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -0.10744 };
	double c2[ 4 ] = { 0.13972, -5.6269e-3, 3.9407e-5, -1.0737e-7 };
	poly_->assign( -9999.9, 76.0, 1, c1 );
	poly_->assign( 76.0, 240.0, 4, c2 );
}
//static double NCPA::N2MoleFractionCalculator::molecular_weight() const  { return 28.0134; }

NCPA::CO2MoleFractionCalculator::CO2MoleFractionCalculator() : NCPA::MoleFractionCalculator(44.0095) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -3.3979 };
	poly_->assign( -9999.9, 240.0, 1, c1 );
}
//static double NCPA::CO2MoleFractionCalculator::molecular_weight() const { return 44.0095; }

NCPA::O3MoleFractionCalculator::O3MoleFractionCalculator() : NCPA::MoleFractionCalculator(47.9982) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -19.027 };
	double c2[ 6 ] = { -19.027, 1.3093, -4.6496e-2, 7.8543e-4, -6.5169e-6, 2.1343e-8 };
	double c3[ 2 ] = { -4.234, -3.0975e-2 };
	poly_->assign( -9999.9, 0.0, 1, c1 );
	poly_->assign( 0.0, 80.0, 6, c2 );
	poly_->assign( 80.0, 240.0, 2, c3 );
}
//static double NCPA::O3MoleFractionCalculator::molecular_weight() const { return 47.9982; }

NCPA::OMoleFractionCalculator::OMoleFractionCalculator() : NCPA::MoleFractionCalculator(15.9994) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -11.195 };
	double c2[ 4 ] = { -11.195, 0.15408, -1.4348e-3, 1.0166e-5 };
	double c3[ 4 ] = { -3.2456, 4.6642e-2, -2.6894e-4, 5.264e-7 };
	poly_->assign( -9999.9, 0.0, 1, c1 );
	poly_->assign( 0.0, 95.0, 4, c2 );
	poly_->assign( 95.0, 240.0, 4, c3 );
}
//static double NCPA::OMoleFractionCalculator::molecular_weight() const { return 15.9994; }

NCPA::NMoleFractionCalculator::NMoleFractionCalculator() : NCPA::MoleFractionCalculator(14.0067) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -53.746 };
	double c2[ 6 ] = { -53.746, 1.5439, -1.8824e-2, 1.1587e-4, -3.5399e-7, 4.2609e-10 };
	poly_->assign( -9999.9, 0.0, 1, c1 );
	poly_->assign( 0.0, 240.0, 6, c2 );
}
//static double NCPA::NMoleFractionCalculator::molecular_weight() const { return 14.0067; }

NCPA::H2OMoleFractionCalculator::H2OMoleFractionCalculator() : NCPA::MoleFractionCalculator(18.0153) {
	poly_ = new NCPA::MultiplePolynomialCalculator<double>();
	double c1[ 1 ] = { -1.7491 };
	double c2[ 6 ] = { -1.7491, 4.4986e-2, -6.8549e-2, 5.4639e-3, -1.5539e-4, 1.5063e-6 };
	double c3[ 6 ] = { -4.2563, 7.6245e-2, -2.1824e-3, -2.3010e-6, 2.4265e-7, -1.2500e-9 };
	double c4[ 6 ] = { -6.2534e-1, -8.3665e-2 };
	poly_->assign( -9999.9, 0.0, 1, c1 );
	poly_->assign( 0.0, 30.0, 6, c2 );
	poly_->assign( 30.0, 100.0, 6, c3 );
	poly_->assign( 100.0, 240.0, 2, c4 );
}
//static double NCPA::H2OMoleFractionCalculator::molecular_weight() const { return 18.0153; }

NCPA::ArMoleFractionCalculator::ArMoleFractionCalculator() : NCPA::MoleFractionCalculator(39.948) {
	n2 = new NCPA::N2MoleFractionCalculator();
}
NCPA::ArMoleFractionCalculator::~ArMoleFractionCalculator() {
	delete n2;
}
double NCPA::ArMoleFractionCalculator::calculate( double z, NCPA::units_t units ) {
	return 0.012 * n2->calculate( z, units );
}
//static double NCPA::ArMoleFractionCalculator::molecular_weight() const { return 39.948; }
