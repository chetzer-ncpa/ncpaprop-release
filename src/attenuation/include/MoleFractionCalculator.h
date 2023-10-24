#ifndef NCPA__ATMOSPHERE_MOLEFRACTIONCALCULATOR_H_INCLUDED
#define NCPA__ATMOSPHERE_MOLEFRACTIONCALCULATOR_H_INCLUDED

#include "NCPACommon.h"
#include "NCPAUnits.h"
#include <utility>
#include <vector>
#include <map>
#include <type_traits>
#include <stdexcept>

namespace NCPA {
	enum class molecular_constituent_t : unsigned int {
		O2,
		N2,
		CO2,
		O3,
		O,
		N,
		H2O,
		AR
	};

	template<typename T>
	class PolynomialCoefficientSet {
	public:
		PolynomialCoefficientSet( std::vector<T> coeffs ) {
			coefficients = coeffs;
		}
		PolynomialCoefficientSet( size_t n, T* coeffs ) {
			coefficients.assign( coeffs, coeffs + n );
		}
		std::vector<T> coefficients;
		T evaluate( T x ) const {
			return NCPA::evalpoly( coefficients, x );
		}
	};

	template<typename T>
	class MultiplePolynomialCalculator {
	public:
		MultiplePolynomialCalculator() {}
		~MultiplePolynomialCalculator() {
			for (auto it = contents.begin(); it != contents.end(); ++it) {
				delete it->second;
			}
		}
		void assign( double range_start, double range_end, size_t n, T *coeffs ) {
			// make sure it's not already there
			if (this->get_set(range_start) != NULL
					|| this->get_set(range_end - std::numeric_limits<double>::epsilon() ) ) {
				throw std::out_of_range("Calculator already contains overlapping range");
			}
			std::pair<double,double> rangepair( range_start, range_end );
			contents[rangepair] = new PolynomialCoefficientSet<T>( n, coeffs );
		}
		T calculate(T value) {
			PolynomialCoefficientSet<T> *set = this->get_set( value );
			if (set != NULL) {
				return set->evaluate( value );
			} else {
				throw std::out_of_range("Coefficient set has not been initialized for this value");
			}
		}
	private:
		std::map< std::pair<double,double>, PolynomialCoefficientSet<T> * > contents;

		PolynomialCoefficientSet<T> * get_set( double x ) {
			for (auto it = contents.begin(); it != contents.end(); ++it) {
				if (x >= it->first.first && x < it->first.second) {
					return it->second;
				}
			}
			return NULL;
		}
	};

	class MoleFractionCalculator {
	public:
		MoleFractionCalculator(double molweight);
		virtual ~MoleFractionCalculator();
		static MoleFractionCalculator * build(molecular_constituent_t type);
		static std::map<NCPA::molecular_constituent_t,NCPA::MoleFractionCalculator*> build_all();
		virtual double calculate( double z, units_t units = NCPA::units_t::DISTANCE_KILOMETERS );
		virtual double molecular_weight() const;
		static std::vector<molecular_constituent_t> available_types();
	protected:
		MultiplePolynomialCalculator<double> *poly_ = NULL;
		const double molecular_weight_;
	};


	class O2MoleFractionCalculator : public MoleFractionCalculator {
	public:
		O2MoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class N2MoleFractionCalculator : public MoleFractionCalculator {
	public:
		N2MoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class CO2MoleFractionCalculator : public MoleFractionCalculator {
	public:
		CO2MoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class O3MoleFractionCalculator : public MoleFractionCalculator {
	public:
		O3MoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class NMoleFractionCalculator : public MoleFractionCalculator {
	public:
		NMoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class OMoleFractionCalculator : public MoleFractionCalculator {
	public:
		OMoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class H2OMoleFractionCalculator : public MoleFractionCalculator {
	public:
		H2OMoleFractionCalculator();
//		virtual static double molecular_weight() const;
	};
	class ArMoleFractionCalculator : public MoleFractionCalculator {
	public:
		ArMoleFractionCalculator();
		~ArMoleFractionCalculator();
		virtual double calculate( double z, units_t units = NCPA::units_t::DISTANCE_KILOMETERS );
//		virtual static double molecular_weight() const { return 39.948; }
	private:
		N2MoleFractionCalculator *n2;
	};

}

#endif
