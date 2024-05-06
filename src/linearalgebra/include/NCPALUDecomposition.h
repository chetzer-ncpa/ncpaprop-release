#ifndef NCPA_LINEARALGEBRA__NCPALUDECOMPOSITION_H_INCLUDED
#define NCPA_LINEARALGEBRA__NCPALUDECOMPOSITION_H_INCLUDED

#include "NCPAMatrix.h"

namespace NCPA {
	template<typename T> class LUDecomposition;
}

/**
 * Swaps two NCPALUDecomposition objects.
 * @param a First object to swap
 * @param b Second object to swap
 */
template<typename T>
void swap( NCPA::LUDecomposition<T> &a, NCPA::LUDecomposition<T> &b ) noexcept;

namespace NCPA {
	template<typename T>
	class LUDecomposition {
	public:
		LUDecomposition() = default;
		LUDecomposition( const NCPA::LUDecomposition<T> &other ) {
			if (other.lower_ != nullptr) {
				this->lower_ = other.lower_->clone();
			}
			if (other.upper_ != nullptr) {
				this->upper_ = other.upper_->clone();
			}
			if (other.permutation_ != nullptr) {
				this->permutation_ = other.permutation_->clone();
			}
		}
		LUDecomposition( NCPA::LUDecomposition<T> &&other ) noexcept {
			::swap( *this, other );
		}

		virtual ~LUDecomposition() {
			this->clear();
		}

		friend void ::swap<T>(
				NCPA::LUDecomposition<T> &a,
				NCPA::LUDecomposition<T> &b ) noexcept;

		virtual NCPA::LUDecomposition<T> * initialize(
				size_t nrows, size_t ncols) = 0;

		virtual NCPA::LUDecomposition<T> * clone() const = 0;

		virtual void compute( const NCPA::Matrix<T> *base,
				bool pivot = true, T tol = ((T)1e-20) ) {

			base->assert_finalized();
			base->assert_square();
			base->assert_sizes_match( this->upper_ );
			base->assert_sizes_match( this->lower_ );
			base->assert_sizes_match( this->permutation_ );

			size_t i, j, k;
			size_t N = base->rows();
			T maxA, absA;

			this->upper_->initialize( N, N );
			this->lower_->initialize( N, N );
			this->permutation_->identity( N, N );
			this->upper_->copy( base );

			for (k = 0; k < N; k++) {
				if (pivot) {
					size_t pivotRow = k;
					T maxVal = this->upper_->get( k, k );

					// find pivot row and swap
					for (i = k+1; i < N; i++) {
						T x = this->upper_->get( i, k );
						if (std::abs(x) > std::abs(maxVal)) {
							maxVal = x;
							pivotRow = i;
						}
					}
					if (std::abs(maxVal) < std::abs(tol)) {
						throw std::runtime_error( "Matrix is degenerate" );
					}
					this->upper_->swap_rows( k, pivotRow );
					this->permutation_->swap_rows( k, pivotRow );
					this->lower_->swap_rows( k, pivotRow );
				}

				// perform elimination
				for (i = k+1; i < N; i++) {
					T ratio = this->upper_->get( i, k ) / this->upper_->get( k, k );
					this->lower_->set( i, k, ratio );
					for (j = 0; j < N; j++) {
						T diff = this->upper_->get( i, j ) - ratio * this->upper_->get( k, j );
						this->upper_->set( i, j, diff );
					}
				}
			}

			for (k = 0; k < N; k++) {
				this->lower_->set( k, k, ((T)1.0) );
			}
			this->lower_->finalize();
			this->upper_->finalize();
			this->permutation_->finalize();
		}
		virtual NCPA::LUDecomposition<T> * clear() {
			if (lower_ != nullptr) {
				delete lower_;
				lower_ = nullptr;
			}
			if (upper_ != nullptr) {
				delete upper_;
				upper_ = nullptr;
			}
			if (permutation_ != nullptr) {
				delete permutation_;
				permutation_ = nullptr;
			}
			return this;
		}


		virtual NCPA::LUDecomposition<T> * initialize(
				const NCPA::Matrix<T> *base ) {
			return this->initialize( base->rows(), base->columns() );
		}

		virtual bool is_ready() const {
			return (lower_ != nullptr)
					&& (upper_ != nullptr)
					&& (permutation_ != nullptr);
		}

		virtual const NCPA::Matrix<T> * lower() const {
			return lower_;
		}
		virtual const NCPA::Matrix<T> * upper() const {
			return upper_;
		}
		virtual const NCPA::Matrix<T> * permutation() const {
			return permutation_;
		}



	protected:
		NCPA::Matrix<T> *lower_ = nullptr,
						*upper_ = nullptr,
						*permutation_ = nullptr;
	};
}

template<typename T>
void swap(
		NCPA::LUDecomposition<T> &a,
		NCPA::LUDecomposition<T> &b ) noexcept {
	using std::swap;
	swap( a.lower_, b.lower_ );
	swap( a.upper_, b.upper_ );
	swap( a.permutation_, b.permutation_ );
}

#endif
