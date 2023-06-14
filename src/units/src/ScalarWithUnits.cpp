#include "units.h"
#include "ScalarWithUnits.h"
#include <utility>

NCPA::ScalarWithUnits::ScalarWithUnits() : value_{ 0.0 }, units_{ NCPA::UNITS_NONE } {}

NCPA::ScalarWithUnits::ScalarWithUnits( double value, units_t units )
	: value_{ value }, units_{ units } {}

NCPA::ScalarWithUnits::ScalarWithUnits( double value, const std::string &units )
	: value_{ value }, units_{ NCPA::Units::fromString( units ) } {}

NCPA::ScalarWithUnits::ScalarWithUnits( const NCPA::ScalarWithUnits &source )
	: value_{ source.value_ }, units_{ source.units_ } {}

NCPA::ScalarWithUnits::ScalarWithUnits( ScalarWithUnits&& that) noexcept {
	swap(*this,that);
}

NCPA::ScalarWithUnits::~ScalarWithUnits() { }

void swap( NCPA::ScalarWithUnits &first, NCPA::ScalarWithUnits &second ) noexcept {
	using std::swap;
	swap(first.value_,second.value_);
	swap(first.units_,second.units_);
}

NCPA::ScalarWithUnits &NCPA::ScalarWithUnits::operator=(NCPA::ScalarWithUnits other) {
	swap(*this, other);
	return *this;
}

NCPA::units_t NCPA::ScalarWithUnits::get_units() const {
	//return units_.top();
	return units_;
}

void NCPA::ScalarWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	//if (new_units != units_.top()) {
	if (new_units != units_) {
		do_units_conversion_( units_, new_units );
	}
	units_ = new_units;
	//units_.push( new_units );
}

void NCPA::ScalarWithUnits::convert_units( const std::string &units ) {
	this->convert_units( NCPA::Units::fromString( units ) );
}

void NCPA::ScalarWithUnits::do_units_conversion_( NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double units_buffer = 0.0;

	// throws out_of_range if conversion is undefined
	units_buffer = NCPA::Units::convert( value_, fromUnits, toUnits );

	// successful, so record the units change
	value_ = units_buffer;
}

double NCPA::ScalarWithUnits::get() const {
	return value_;
}

double NCPA::ScalarWithUnits::get_as( NCPA::units_t u ) const {
	return NCPA::Units::convert( value_, units_, u );
}

double NCPA::ScalarWithUnits::get_as( const std::string &u ) const {
	return this->get_as( NCPA::Units::fromString( u ) );
}

void NCPA::ScalarWithUnits::set_value( double newval ) {
	value_ = newval;
}

void NCPA::ScalarWithUnits::set_units( NCPA::units_t new_units ) {
	units_ = new_units;
}

void NCPA::ScalarWithUnits::set_units( const std::string &new_units ) {
	this->set_units( NCPA::Units::fromString( new_units ) );
}

void NCPA::ScalarWithUnits::set( double newval, NCPA::units_t new_units ) {
	set_value( newval );
	set_units( new_units );
}

void NCPA::ScalarWithUnits::set( double newval, const std::string &new_units ) {
	this->set( newval, NCPA::Units::fromString( new_units ) );
}


std::ostream &NCPA::operator<<( std::ostream &output, const NCPA::ScalarWithUnits &D ) {
	output << D.get() << " " << NCPA::Units::toStr( D.get_units() );
	return output;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator+(
		const NCPA::ScalarWithUnits &second ) const {
	NCPA::ScalarWithUnits temp1(*this);
	temp1 += second;
	return temp1;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator-(
		const NCPA::ScalarWithUnits &second ) const {
	NCPA::ScalarWithUnits temp1(*this);
	temp1 -= second;
	return temp1;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator+=( NCPA::ScalarWithUnits const &second ) {
	NCPA::ScalarWithUnits temp2(second);
	temp2.convert_units( get_units() );
	value_ += temp2.get();
	return *this;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator-=( NCPA::ScalarWithUnits const &second ) {
	NCPA::ScalarWithUnits temp2(second);
	temp2.convert_units( get_units() );
	value_ -= temp2.get();
	return *this;
}

