#ifndef NCPA_VECTORWITHUNITS_H_INCLUDED
#define NCPA_VECTORWITHUNITS_H_INCLUDED

#include "units.h"
#include "ScalarWithUnits.h"
#include <vector>

namespace NCPA { class VectorWithUnits; }
void swap( NCPA::VectorWithUnits&, NCPA::VectorWithUnits& ) noexcept;


namespace NCPA {
	class VectorWithUnits : public std::vector<ScalarWithUnits> {
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

			// methods
			virtual void as_array( NCPA::ScalarWithUnits *&buffer, bool normFirst=true );
			virtual void as_array( double *&buffer, units_t &units, bool normFirst=true );
			virtual void as_array( double *&buffer, bool normFirst=true );

			virtual void convert_units( units_t new_units );
			virtual void convert_units( const std::string &new_units );

			// Fill the vector with identical values.  Does not resize.
			virtual void fill( double value, units_t units );
			virtual void fill( double value, const std::string &units );
			virtual void fill( const ScalarWithUnits &value );

			virtual units_t get_units( bool normFirst=true );
			virtual void get_values( size_t &n, double* buffer, bool normFirst=true );
			virtual void get_values( double* buffer, bool normFirst=true );

			virtual bool is_normalized() const;
			virtual void normalize_units();

			virtual void set( size_t n_points, const double *values, units_t units );
			virtual void set( size_t n_points, const double *values, const std::string &units );
			virtual void set( size_t n_points, const ScalarWithUnits *values );

			virtual void set_units( units_t new_units );
			virtual void set_units( const std::string &new_units );
	};

}

#endif
