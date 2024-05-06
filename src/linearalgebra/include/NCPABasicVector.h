#ifndef NCPA_LINEARALGEBRA__NCPABASICVECTOR_H_INCLUDED
#define NCPA_LINEARALGEBRA__NCPABASICVECTOR_H_INCLUDED

#include "NCPACommon.h"
#include <map>
#include <vector>
#include <stdexcept>
#include <sstream>

namespace NCPA { template<typename T> class BasicVector; }

template<typename T>
void swap( NCPA::BasicVector<T> &a, NCPA::BasicVector<T> &b ) noexcept;

namespace NCPA {

	/**
	 * BasicVector class.
	 * Provides an interface similar to that of the std::vector class, but with
	 * underlying storage as an map to optimize for situations where
	 * few of the vector values are nonzero.
	 */
	template<typename T>
	class BasicVector : public std::map<size_t,T>, public NCPA::Vector<T> {
	public:

		/**
		 * Default constructor.
		 * Generates an empty vector.
		 */
		BasicVector() : std::map<size_t,T>(), NCPA::Vector<T>(),
		vector_length_{0} {}

		/**
		 * Allocation constructor.
		 * Creates an empty vector of size n.
		 * @param n The size of the vector (not the number of nonzero elements)
		 */
		BasicVector( size_t n ) : std::map<size_t,T>(), NCPA::Vector<T>() ,
				vector_length_{n} {}

		/**
		 * Initialization constructor.
		 * Assumes the values in <code>vals</code> are complete and in order.
		 * @param vals The values.
		 */
		BasicVector( std::vector<T> &vals ) : std::map<size_t,T>(),
				NCPA::Vector<T>(), vector_length_{vals.size()} {

//			this->reserve( vals.size() );
			for (size_t i = 0; i < vals.size(); i++) {
				(*this)[i] = vals[i];
			}
		}

		/**
		 * Initialization constructor.
		 * Associates the indices in <code>inds</code> with the values in
		 * <code>vals</code>.
		 * @param n The size of the vector
		 * @param inds The indices to store.
		 * @param vals The values associated with the indices.
		 * @throws std::invalid_argument if <code>inds</code> and <code>vals</code> are
		 * 		different sizes.
		 */
		BasicVector( size_t n, std::vector<size_t> inds, std::vector<T> vals ) :
			std::map<size_t,T>(), NCPA::Vector<T>(), vector_length_{n} {

//			this->reserve( vals.size() );
			for (size_t i = 0; i < inds.size(); i++) {
				(*this)[ inds[i] ] = vals[ i ];
			}
		}

		/**
		 * Initialization constructor.
		 * Associates the indices in <code>inds</code> with the values in
		 * <code>vals</code>.
		 * @param n The size of the vector.
		 * @param nvals The number of values to store.
		 * @param inds The indices to store.
		 * @param vals The values associated with the indices.
		 */
		BasicVector( size_t n, size_t nvals, const size_t *inds, const T *vals ) :
			NCPA::BasicVector<T>(
					n,
					std::vector<size_t>( inds, inds+nvals ),
					std::vector<T>( vals, vals+nvals ) )
			{}

		/**
		 * Initialization constructor.
		 * Associates the indices in <code>inds</code> with the values in
		 * <code>vals</code>.
		 * @param n The size of the vector.
		 * @param inds The indices to store.
		 * @param vals The values associated with the indices.
		 * @throws std::invalid_argument if <code>inds</code> and <code>vals</code> are
		 * 		different sizes.
		 */
		BasicVector( size_t n, std::initializer_list<size_t> inds, std::initializer_list<T> vals) :
			NCPA::BasicVector<T>( n, std::vector<size_t>( inds ), std::vector<T>( vals ) )
			{}

		/**
		 * Constant initialization constructor.
		 * Associates all indices with the value <code>val</code>.
		 * @param n The size of the vector
		 * @param val The value to store at each index.
		 */
		BasicVector( size_t n, T val ) : NCPA::BasicVector<T>( n ) {
			for (size_t i = 0; i < n; i++) {
				this->set( i, val );
			}
		}


		/**
		 * Constant initialization constructor.
		 * Associates the indices in <code>inds</code> with the value <code>val</code>.
		 * @param n The size of the vector
		 * @param inds The indices to store.
		 * @param val The value to store at each index.
		 */
		BasicVector( size_t n, std::vector<size_t> inds, T val ) :
			NCPA::BasicVector<T>( n, inds, std::vector<T>( inds.size(), val ) )
			{}

		/**
		 * Constant initialization constructor.
		 * Associates the indices in <code>inds</code> with the value <code>val</code>.
		 * @param n The size of the vector
		 * @param nvals The number of values to store.
		 * @param inds The indices to store.
		 * @param val The value to store at each index.
		 */
		BasicVector( size_t n, size_t nvals, const size_t *inds, T val ) :
			NCPA::BasicVector<T>(
				n,
				std::vector<size_t>( inds, inds+nvals ),
				std::vector<T>( nvals, val ) )
			{}

		/**
		 * Constant initialization constructor.
		 * Associates the indices in <code>inds</code> with the value <code>val</code>.
		 * @param n The size of the vector
		 * @param inds The indices to store.
		 * @param val The value to store at each index.
		 */
		BasicVector( size_t n, std::initializer_list<size_t> inds, T val ) :
			NCPA::BasicVector<T>( n, std::vector<size_t>( inds ),
					std::vector<T>( inds.size(), val ) )
			{}


		/**
		 * Copy constructor.
		 * @param source The vector to copy.
		 */
		BasicVector( const NCPA::BasicVector<T> &source ) :
			std::map<size_t,T>( source ), NCPA::Vector<T>( source ) {
			this->vector_length_ = source.vector_length_;
		}
//
//		/**
//		 * Superclass copy constructor.
//		 * @param source The vector to copy
//		 */
//		BasicVector( const NCPA::Vector<T> *source ) :
//			std::map<size_t,T>(), NCPA::Vector<T>( *source ),
//			vector_length_{ source->size() } {
//
//		}

		/**
		 * Move constructor.
		 * @param source The vector to assimilate.
		 */
		BasicVector( NCPA::BasicVector<T> &&source ) noexcept :
			NCPA::BasicVector<T>() {
			::swap( *this, source );
		}

		virtual ~BasicVector() {
			this->clear();
		}

		/**
		 * Assignment operator.
		 * @param other The vector to assign to this.
		 */
		BasicVector<T>& operator=( NCPA::BasicVector<T> other ) {
			::swap(*this,other);
			return *this;
		}

		friend void ::swap<T>( NCPA::BasicVector<T> &a, NCPA::BasicVector<T> &b ) noexcept;

		virtual NCPA::Vector<T> * clear() {
			static_cast<std::map<size_t,T>*>( this )->clear();
			this->vector_length_ = 0;
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual NCPA::Vector<T> * clear( size_t ind ) {
			auto it = this->find( ind );
			if (it != this->end()) {
				this->erase( it );
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual NCPA::Vector<T> * resize( size_t n ) {
			if (n < vector_length_) {
				NCPA::BasicVector<T> *temp = new NCPA::BasicVector<T>( n );
				for (auto it = this->begin(); it != this->end(); ++it) {
					if (it->first < n) {
						temp->set( it->first, it->second );
					}
				}
				swap<T>( *this, *temp );
				delete temp;
			} else {
				vector_length_ = n;
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual size_t size() const {
			return vector_length_;
		}

		virtual NCPA::Vector<T> * clone() const {
			return static_cast<NCPA::Vector<T> *>( new NCPA::BasicVector<T>( *this ) );
		}

		virtual T get( size_t n ) const {
			this->assert_element_in_range( n );
			auto it = this->find( n );
			if (it == this->cend()) {
				return ((T)0);
			} else {
				return it->second;
			}
		}

		virtual void as_array( size_t n, T *&vals ) const {
			this->assert_size_matches( n );
			if (vals == nullptr) {
				vals = NCPA::zeros<T>( n );
			}
			for (auto it = this->cbegin(); it != this->cend(); ++it) {
				vals[it->first] = it->second;
			}
		}

		virtual NCPA::Vector<T> * set( size_t n, T val ) {
			this->assert_element_in_range( n );
			(*this)[ n ] = val;
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual NCPA::Vector<T> * scale( T val ) {
			for (auto it = this->begin(); it != this->end(); ++it) {
				it->second *= val;
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual NCPA::Vector<T> * add( const NCPA::Vector<T> *b ) {
			this->assert_sizes_match( *b );
			auto binds = b->get_defined_indices();
			for (auto it = binds.cbegin(); it != binds.cend(); ++it) {
				this->set( *it, this->get(*it) + b->get(*it) );
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual NCPA::Vector<T> * subtract( const NCPA::Vector<T> *b ) {
			this->assert_sizes_match( *b );
			auto binds = b->get_defined_indices();
			for (auto it = binds.cbegin(); it != binds.cend(); ++it) {
				this->set( *it, this->get(*it) - b->get(*it) );
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual T scalar_product( const NCPA::Vector<T> *b ) const {
			this->assert_sizes_match( *b );
			T product = ((T)0);
			for (auto it = this->cbegin(); it != this->cend(); ++it) {
				product += it->second * b->get( it->first );
			}
			return product;
		}

		virtual NCPA::Vector<T> * product( const NCPA::Vector<T> *b ) const {
			this->assert_sizes_match( *b );
			NCPA::Vector<T> *prod = this->clone();
			for (auto it = this->cbegin(); it != this->cend(); ++it) {
				prod->set( it->first, it->second * b->get( it->first ) );
			}
			return prod;
		}

		virtual std::vector<T> as_std() const {
			std::vector<T> v( vector_length_ );
			for (auto it = this->cbegin(); it != this->cend(); ++it) {
				v[ it->first ] = it->second;
			}
			return v;
		}

		virtual NCPA::Vector<T> * from_std( const std::vector<T> &from ) {
			this->clear();
//			this->reserve( from.size() );
			this->vector_length_ = from.size();
			for (size_t i = 0; i < from.size(); i++) {
				this->set( i, from[i] );
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual NCPA::Vector<T> * from_array( size_t n, const T *vals ) {
			this->clear();
			this->vector_length_ = n;
			for (size_t i = 0; i < n; i++) {
				this->set( i, vals[i] );
			}
			return static_cast<NCPA::Vector<T> *>( this );
		}

		virtual std::vector<size_t> get_defined_indices() const {
			std::vector<size_t> inds;
			inds.reserve( this->size() );
			for (auto it = this->cbegin(); it != this->cend(); ++it) {
				inds.push_back( it->first );
			}
//			std::sort( inds.begin(), inds.end() );
			return inds;
		}

		virtual void copy( const NCPA::Vector<T> *other ) {
			this->clear();
			this->resize( other->size() );
			std::vector<size_t> inds = other->get_defined_indices();
			for (auto it = inds.cbegin(); it != inds.cend(); it++) {
				this->set( *it, other->get( *it ) );
			}
		}

	protected:
		size_t vector_length_ = 0;
	};

}


template<typename T>
void swap( NCPA::BasicVector<T> &a, NCPA::BasicVector<T> &b ) noexcept {
	using std::swap;
//	::swap<T>( static_cast<NCPA::Vector<T> &>(a), static_cast<NCPA::Vector<T> &>(b) );
	swap( static_cast<std::map<size_t,T> &>(a),
			static_cast<std::map<size_t,T> &>(b) );
	swap( a.vector_length_, b.vector_length_ );
}

#endif




