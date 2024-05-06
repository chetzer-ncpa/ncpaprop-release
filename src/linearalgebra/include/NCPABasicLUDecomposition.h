#ifndef NCPA_LINEARALGEBRA__NCPABASICLUDECOMPOSITION_H_INCLUDED
#define NCPA_LINEARALGEBRA__NCPABASICLUDECOMPOSITION_H_INCLUDED

#include "NCPAMatrix.h"
#include "NCPABasicMatrix.h"
#include "NCPALUDecomposition.h"

namespace NCPA {
	template<typename T> class BasicLUDecomposition;
}

/**
 * Swaps two NCPABasicLUDecomposition objects.
 * @param a First object to swap
 * @param b Second object to swap
 */
template<typename T>
void swap( NCPA::BasicLUDecomposition<T> &a, NCPA::BasicLUDecomposition<T> &b ) noexcept;

namespace NCPA {
	template<typename T>
	class BasicLUDecomposition : public NCPA::LUDecomposition<T> {
	public:
		BasicLUDecomposition() = default;
		BasicLUDecomposition(
				const NCPA::BasicLUDecomposition<T> &other )
				: NCPA::LUDecomposition<T>( other ) {}
		BasicLUDecomposition( NCPA::BasicLUDecomposition<T> &&other )
				noexcept {
			::swap( *this, other );
		}
		virtual ~BasicLUDecomposition() {
			this->clear();
		}

		friend void ::swap<T>(
				NCPA::BasicLUDecomposition<T> &a,
				NCPA::BasicLUDecomposition<T> &b ) noexcept;
		BasicLUDecomposition<T>& operator=(
				NCPA::BasicLUDecomposition<T> other ) {
			::swap(*this,other);
			return *this;
		}

		virtual NCPA::LUDecomposition<T> * clone() const {
			return static_cast<NCPA::LUDecomposition<T> *>(
					new NCPA::BasicLUDecomposition<T>( *this ) );
		}

		virtual NCPA::LUDecomposition<T> * initialize(
				size_t nrows, size_t ncols) {
			this->clear();
			this->lower_ = new NCPA::BasicMatrix<T>( nrows, ncols );
			this->upper_ = new NCPA::BasicMatrix<T>( nrows, ncols );
			this->permutation_ = new NCPA::BasicMatrix<T>( nrows, ncols );
			return static_cast<NCPA::LUDecomposition<T> *>( this );
		}
	};
}

template<typename T>
void swap(
		NCPA::BasicLUDecomposition<T> &a,
		NCPA::BasicLUDecomposition<T> &b ) noexcept {
	using std::swap;
	::swap(
			static_cast<NCPA::LUDecomposition<T> &>(a),
			static_cast<NCPA::LUDecomposition<T> &>(b)
			);
}

#endif
