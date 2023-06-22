#include "units.h"
#include "ScalarWithUnits.h"
#include <utility>
#include <stdexcept>

NCPA::ScalarWithUnits::ScalarWithUnits() : value_{ 0.0 }, units_{ NCPA::units_t::NONE } {}

NCPA::ScalarWithUnits::ScalarWithUnits( double value )
	: value_{ value }, units_{ NCPA::units_t::NONE } {}

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

std::ostream &operator<<( std::ostream &output, const NCPA::ScalarWithUnits &D ) {
	output << D.get() << " " << NCPA::Units::toStr( D.get_units() );
	return output;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator+=( NCPA::ScalarWithUnits const &second ) {
	if (NCPA::Units::can_convert(second.units_, this->units_)) {
		this->value_ += second.get_as( this->units_ );
		return *this;
	} else {
		throw NCPA::invalid_conversion(this->units_,second.units_);
	}
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

//NCPA::ScalarWithUnits operator+( double first, const NCPA::ScalarWithUnits &second ) {
//	return NCPA::ScalarWithUnits(second+first);
//}
//
//NCPA::ScalarWithUnits operator-( double first, const NCPA::ScalarWithUnits &second ) {
//	return NCPA::ScalarWithUnits( (-second) + first );
//}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator-=( NCPA::ScalarWithUnits const &second ) {
	*this += -second;
	return *this;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator*=( double second ) {
	this->value_ *= second;
	return *this;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator/=( double second ) {
	this->value_ /= second;
	return *this;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator*( double second ) {
	NCPA::ScalarWithUnits temp(*this);
	temp *= second;
	return temp;
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator/( double second ) {
	NCPA::ScalarWithUnits temp(*this);
	temp /= second;
	return temp;
}

double NCPA::ScalarWithUnits::over( const NCPA::ScalarWithUnits &second ) const {
	return this->value_ / second.as(this->units_);
}

NCPA::ScalarWithUnits NCPA::ScalarWithUnits::operator-() const {
	NCPA::ScalarWithUnits temp1(*this);
	temp1.value_ = -temp1.value_;
	return temp1;
}

bool operator==( const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b ) {
	if (NCPA::Units::can_convert(a.units_, b.units_)) {
		NCPA::ScalarWithUnits bb = b;
		bb.convert_units( a.units_ );
		return (a.value_ == b.value_);
	} else {
		return false;
	}
}

bool operator!=( const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b ) {
	return !(a == b);
}

bool operator>( const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b ) {
	if (NCPA::Units::can_convert(b.units_, a.units_)) {
		NCPA::ScalarWithUnits bb = b;
		bb.convert_units( a.units_ );
		return (a.value_ > b.value_);
	} else {
		throw NCPA::invalid_conversion(a.units_,b.units_);
	}
}

bool operator<( const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b ) {
	if (NCPA::Units::can_convert(b.units_, a.units_)) {
		NCPA::ScalarWithUnits bb = b;
		bb.convert_units( a.units_ );
		return (a.value_ < b.value_);
	} else {
		throw NCPA::invalid_conversion(a.units_,b.units_);
	}
}

bool operator>=( const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b ) {
	return !(a < b);
}

bool operator<=( const NCPA::ScalarWithUnits &a, const NCPA::ScalarWithUnits &b ) {
	return !(a > b);
}




double NCPA::ScalarWithUnits::as( NCPA::units_t u ) const {
	return this->get_as(u);
}

double NCPA::ScalarWithUnits::as( const std::string &u ) const {
	return this->get_as(u);
}

void NCPA::ScalarWithUnits::convert( NCPA::units_t new_units ) {
	this->convert_units( new_units );
}

void NCPA::ScalarWithUnits::convert( const std::string &new_units ) {
	this->convert_units( new_units );
}

void NCPA::ScalarWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw invalid_conversion and leave original units unchanged if there's an error
	if (new_units != units_) {
		do_units_conversion_( units_, new_units );
	}
	units_ = new_units;
}

void NCPA::ScalarWithUnits::convert_units( const std::string &units ) {
	this->convert_units( NCPA::Units::fromString( units ) );
}

void NCPA::ScalarWithUnits::do_units_conversion_( NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double units_buffer = 0.0;

	// throws invalid_conversion if conversion is undefined
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






















NCPA::units_t NCPA::ScalarWithUnits::get_units() const {
	//return units_.top();
	return units_;
}




void NCPA::ScalarWithUnits::set( double newval, NCPA::units_t new_units ) {
	set_value( newval );
	set_units( new_units );
}

void NCPA::ScalarWithUnits::set( double newval, const std::string &new_units ) {
	this->set( newval, NCPA::Units::fromString( new_units ) );
}

void NCPA::ScalarWithUnits::set_units( NCPA::units_t new_units ) {
	units_ = new_units;
}

void NCPA::ScalarWithUnits::set_units( const std::string &new_units ) {
	this->set_units( NCPA::Units::fromString( new_units ) );
}

void NCPA::ScalarWithUnits::set_value( double newval ) {
	value_ = newval;
}


