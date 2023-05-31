#ifndef NCPA_VECTORWITHUNITS_H_INCLUDED
#define NCPA_VECTORWITHUNITS_H_INCLUDED

#include "units.h"
#include "ScalarWithUnits.h"


namespace NCPA {
	class VectorWithUnits {
		protected:
			size_t n_;
			double *values_;
			units_t units_;

			void do_units_conversion_( size_t n_points, double *inplace,
				NCPA::units_t fromUnits, NCPA::units_t toUnits );

		public:
			VectorWithUnits();
			VectorWithUnits( size_t n_points, const double *values, units_t units );
			VectorWithUnits( size_t n_points, const double *values, const std::string &units );
			VectorWithUnits( size_t n_points, const ScalarWithUnits *values );
			VectorWithUnits( const VectorWithUnits &source );
			virtual ~VectorWithUnits();

			virtual void convert_units( units_t new_units );
			virtual void convert_units( const std::string &new_units );
			virtual units_t get_units() const;
			//virtual void revert_units();

			virtual size_t size() const;
			virtual void get_vector( double *buffer, units_t *buffer_units ) const;
			virtual void get_vector( double *buffer ) const;

			// operator overloading
			double operator []( size_t i ) const;
			double & operator []( size_t i );
		};

}

#endif
