#ifndef NCPA_VECTORWITHUNITS_H_INCLUDED
#define NCPA_VECTORWITHUNITS_H_INCLUDED

#include "units.h"
#include "ScalarWithUnits.h"
#include <vector>

namespace NCPA { class VectorWithUnits; }
void swap( NCPA::VectorWithUnits&, NCPA::VectorWithUnits& ) noexcept;


namespace NCPA {
	class VectorWithUnits : public std::vector<double> {
		protected:
			units_t units_;

		public:
			// constructors
			VectorWithUnits();
			VectorWithUnits( size_t n_points, const double *values, units_t units );
			VectorWithUnits( size_t n_points, const double *values, const std::string &units );
			VectorWithUnits( size_t n_points, const ScalarWithUnits *values );
			VectorWithUnits( size_t n_points, const ScalarWithUnits &singleValue );
			VectorWithUnits( size_t n_points, double singleValue, units_t units );

			// copy constructor
			VectorWithUnits( const VectorWithUnits &source );

			// move constructor
			VectorWithUnits( VectorWithUnits &&source ) noexcept;

			// destructor
			virtual ~VectorWithUnits();

			// assignment and swapping
			friend void ::swap( VectorWithUnits &first, VectorWithUnits &second ) noexcept;
			VectorWithUnits &operator=(VectorWithUnits other);

			virtual void convert_units( units_t new_units );
			virtual void convert_units( const std::string &new_units );
			virtual units_t get_units() const;
			virtual void set( size_t n_points, const double *values, NCPA::units_t units );
			virtual void set_units( units_t new_units );
			virtual void fill( size_t n_points, double value );
			virtual void fill( double value );

			virtual void set_values( size_t n_points, const double *values );
			virtual void get_vector( double *buffer, units_t *buffer_units ) const;
			virtual void get_vector( double *buffer ) const;


		};

}

#endif
