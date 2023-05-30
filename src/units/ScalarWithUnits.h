#ifndef NCPA_SCALARWITHUNITS_H_INCLUDED
#define NCPA_SCALARWITHUNITS_H_INCLUDED

#include "units.h"



namespace NCPA {
	class ScalarWithUnits {
		protected:
			double value_;
			//std::stack< NCPA::units_t > units_;
			units_t units_;
			void do_units_conversion_( NCPA::units_t fromUnits, NCPA::units_t toUnits );

		public:
			ScalarWithUnits();
			ScalarWithUnits( double value, units_t property_units );
			ScalarWithUnits( double value, const std::string &units );
			ScalarWithUnits( const ScalarWithUnits &source );
			virtual ~ScalarWithUnits();

			virtual double get() const;
			virtual double get_as( units_t u ) const;
			virtual double get_as( const std::string &units ) const;
			virtual void set_value( double newval );
			virtual void set_units( units_t new_units );
			virtual void set_units( const std::string &units );
			virtual void set( double newval, units_t new_units );
			virtual void set( double newval, const std::string &units );

			virtual void convert_units( units_t new_units );
			virtual void convert_units( const std::string &units );
			virtual units_t get_units() const;
			//virtual void revert_units();

			friend ScalarWithUnits operator+( ScalarWithUnits first, ScalarWithUnits const &second );
			friend ScalarWithUnits operator-( ScalarWithUnits first, ScalarWithUnits const &D );
			ScalarWithUnits operator+=( ScalarWithUnits const &second );
			ScalarWithUnits operator-=( ScalarWithUnits const &second );


		};
		std::ostream &operator<<( std::ostream &output, const ScalarWithUnits &D );

}

#endif
