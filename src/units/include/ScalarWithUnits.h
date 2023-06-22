#ifndef NCPA_SCALARWITHUNITS_H_INCLUDED
#define NCPA_SCALARWITHUNITS_H_INCLUDED

#include "units.h"

// forward declarations
namespace NCPA { class ScalarWithUnits; }
void swap( NCPA::ScalarWithUnits&, NCPA::ScalarWithUnits& ) noexcept;
bool operator==(const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b);
bool operator!=(const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b);
bool operator>=(const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b);
bool operator<=(const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b);
bool operator>(const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b);
bool operator<(const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b);

std::ostream &operator<<( std::ostream &output, const NCPA::ScalarWithUnits &D );

//bool operator==(const NCPA::ScalarWithUnits &a, double b);

// class definitions
namespace NCPA {
	class ScalarWithUnits {
		public:
			ScalarWithUnits();
			ScalarWithUnits( double value );
			ScalarWithUnits( double value, units_t property_units );
			ScalarWithUnits( double value, const std::string &units );

			// copy constructor
			ScalarWithUnits( const ScalarWithUnits &source );

			// move constructor
			ScalarWithUnits( ScalarWithUnits&& that) noexcept;

			// destructor
			virtual ~ScalarWithUnits();

			// swap
			friend void ::swap(ScalarWithUnits &first, ScalarWithUnits &second) noexcept;

			virtual double get() const;
			virtual units_t get_units() const;
			virtual double get_as( units_t u ) const;
			virtual double get_as( const std::string &units ) const;
			virtual double as( units_t u ) const;
			virtual double as( const std::string &units ) const;

			virtual void set_value( double newval );
			virtual void set_units( units_t new_units );
			virtual void set_units( const std::string &units );
			virtual void set( double newval, units_t new_units );
			virtual void set( double newval, const std::string &units );

			virtual void convert( units_t new_units );
			virtual void convert( const std::string &units );
			virtual void convert_units( units_t new_units );
			virtual void convert_units( const std::string &units );

			ScalarWithUnits &operator=(ScalarWithUnits other);
			ScalarWithUnits operator+( const ScalarWithUnits &second ) const;
			ScalarWithUnits operator-( const ScalarWithUnits &second ) const;
			ScalarWithUnits operator+=( ScalarWithUnits const &second );
			ScalarWithUnits operator-=( ScalarWithUnits const &second );
			ScalarWithUnits operator*=( double factor );
			ScalarWithUnits operator/=( double factor );
			ScalarWithUnits operator*( double factor );
			ScalarWithUnits operator/( double factor );
			ScalarWithUnits operator-() const;

			double over( const ScalarWithUnits &b ) const;

			friend bool ::operator==(const ScalarWithUnits &a, const ScalarWithUnits &b);
			friend bool ::operator!=(const ScalarWithUnits &a, const ScalarWithUnits &b);
			friend bool ::operator>=(const ScalarWithUnits &a, const ScalarWithUnits &b);
			friend bool ::operator<=(const ScalarWithUnits &a, const ScalarWithUnits &b);
			friend bool ::operator>(const ScalarWithUnits &a, const ScalarWithUnits &b);
			friend bool ::operator<(const ScalarWithUnits &a, const ScalarWithUnits &b);
			friend std::ostream &operator<<( std::ostream &output, const ScalarWithUnits &D );

		protected:
			double value_;
			//std::stack< NCPA::units_t > units_;
			units_t units_;
			void do_units_conversion_( NCPA::units_t fromUnits, NCPA::units_t toUnits );



		};

//		ScalarWithUnits operator+( ScalarWithUnits first, ScalarWithUnits const &second );
//		ScalarWithUnits operator-( ScalarWithUnits first, ScalarWithUnits const &D );

}

#endif
