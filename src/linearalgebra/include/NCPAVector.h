#ifndef NCPA__LINEARALGEBRA_NCPAVECTOR_H_INCLUDED_
#define NCPA__LINEARALGEBRA_NCPAVECTOR_H_INCLUDED_

#include <vector>
#include <sstream>
#include <stdexcept>

namespace NCPA {
	template<typename T> class Vector;
}

/**
 * Swaps two Vector objects.
 * @param a First object to swap
 * @param b Second object to swap
 */
template<typename T>
void swap( NCPA::Vector<T> &a, NCPA::Vector<T> &b ) noexcept;

template<typename T>
bool operator==(const NCPA::Vector<T> &a, const NCPA::Vector<T> &b);
template<typename T>
bool operator!=(const NCPA::Vector<T> &a, const NCPA::Vector<T> &b);

namespace NCPA {

//	/**
//	 * An enumerated type for different kinds of matrices.
//	 */
//	enum class vector_t : unsigned int {
//		UNDEFINED = 0,		/**< Generic, undefined vector */
//		DENSE,				/**< Standard dense (i.e. not storage-optimized) vector */
//		SPARSE				/**< Sparse, storage-optimized vector */
//	};

	template<typename T>
	class Vector {

	// prevent size_t and the like, matrix elements must be capable of being signed
	static_assert(!(std::is_unsigned<T>::value),
			"Vector cannot be instantiated with unsigned arithmetic types.");

	public:
		Vector() {}
		Vector( const NCPA::Vector<T> &other ) {}
		Vector( NCPA::Vector<T> &&other ) noexcept {
			::swap( *this, other );
		}
		virtual ~Vector() {}
		friend void ::swap<T>( NCPA::Vector<T> &a, NCPA::Vector<T> &b ) noexcept;

		/**
		 * @defgroup ncpa-vector-api NCPA::Vector API
		 * @{
		 */
		virtual NCPA::Vector<T> * clear() = 0;
		virtual NCPA::Vector<T> * clear( size_t ind ) = 0;

		virtual NCPA::Vector<T> * resize( size_t n ) = 0;

		virtual NCPA::Vector<T> * clone() const = 0;

		virtual size_t size() const = 0;

		virtual T get( size_t n ) const = 0;

		virtual void as_array( size_t n, T *&vals ) const = 0;

		virtual NCPA::Vector<T> * set( size_t n, T val ) = 0;

		virtual NCPA::Vector<T> * scale( T val ) = 0;

		virtual NCPA::Vector<T> * add( const NCPA::Vector<T> *b ) = 0;

		virtual NCPA::Vector<T> * subtract( const NCPA::Vector<T> *b ) = 0;

		virtual T scalar_product( const NCPA::Vector<T> *b ) const = 0;

		virtual NCPA::Vector<T> * product( const NCPA::Vector<T> *b ) const = 0;

		virtual std::vector<T> as_std() const = 0;

		virtual NCPA::Vector<T> * from_std( const std::vector<T> &v ) = 0;

		virtual NCPA::Vector<T> * from_array( size_t n, const T *vals ) = 0;

		virtual std::vector<size_t> get_defined_indices() const = 0;

		virtual void copy( const NCPA::Vector<T> *other ) = 0;

		/**
		 * @}
		 */

		virtual void assert_size_matches( size_t n ) const {
			if (this->size() != n) {
				std::ostringstream oss;
				oss 	<< "Vector expected to be size " << n << " but is of size "
						<< this->size();
				throw std::out_of_range( oss.str() );
			}
		}

		virtual void assert_sizes_match( const NCPA::Vector<T> &b ) const {
			if (this->size() != b.size()) {
				std::ostringstream oss;
				oss 	<< "Vector expected to be size " << b.size() << " but is of size "
						<< this->size();
				throw std::out_of_range( oss.str() );
			}
		}

		virtual void assert_element_in_range( size_t n ) const {
			if (this->size() <= n) {
				std::ostringstream oss;
				oss 	<< "Element " << n << " is out of range for vector of size "
						<< this->size();
				throw std::out_of_range( oss.str() );
			}
		}

		virtual bool is_unequal( const NCPA::Vector<T> *other ) const {
			if (this->size() != other->size()) {
				return true;
			}
			for (size_t i = 0; i < this->size(); i++) {
				if (this->get(i) != other->get(i)) {
					return true;
				}
			}
			return false;
		}
	};
}


template<typename T>
void swap( NCPA::Vector<T> &a, NCPA::Vector<T> &b ) noexcept {}


#endif
