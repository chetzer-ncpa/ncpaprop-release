#ifndef NCPA_LINEARALGEBRA__NCPABASICLINEARALGEBRAFACTORY_H_INCLUDED
#define NCPA_LINEARALGEBRA__NCPABASICLINEARALGEBRAFACTORY_H_INCLUDED

#include <vector>

#include "NCPALinearAlgebraFactory.h"

// vectors
#include "NCPAVector.h"
#include "NCPABasicVector.h"

// matrices
#include "NCPAMatrix.h"
#include "NCPABasicMatrix.h"

// Decompositions
#include "NCPALUDecomposition.h"
#include "NCPABasicLUDecomposition.h"

// Solvers



namespace NCPA {
	template<typename T>
	class BasicLinearAlgebraFactory : public NCPA::LinearAlgebraFactory<T> {
	public:
		BasicLinearAlgebraFactory() = default;
		virtual ~BasicLinearAlgebraFactory() {}

		// vectors
		virtual NCPA::Vector<T>* build_vector() const {
			return new NCPA::BasicVector<T>();
		}
		virtual NCPA::Vector<T>* build_vector( size_t nvals ) const {
			return new NCPA::BasicVector<T>( nvals );
		}
		virtual NCPA::Vector<T>* build_vector( std::vector<T> vals ) const {
			return new NCPA::BasicVector<T>( vals );
		}
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::vector<size_t> inds, std::vector<T> vals ) const {
			return new NCPA::BasicVector<T>( n, inds, vals );
		}
		virtual NCPA::Vector<T>* build_vector( size_t n,
				size_t nvals, const size_t *inds, const T *vals ) const {
			return new NCPA::BasicVector<T>( n, nvals, inds, vals );
		}
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::initializer_list<size_t> inds,
				std::initializer_list<T> vals ) const {
			return new NCPA::BasicVector<T>( n, inds, vals );
		}
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::vector<size_t> inds, T val ) const {
			return new NCPA::BasicVector<T>( n, inds, val );
		}
		virtual NCPA::Vector<T>* build_vector( size_t n, size_t nvals,
				const size_t *inds, T val ) const {
			return new NCPA::BasicVector<T>( n, nvals, inds, val );
		}
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::initializer_list<size_t> inds, T val ) const {
			return new NCPA::BasicVector<T>( n, inds, val );
		}

		// matrices
		virtual NCPA::Matrix<T>* build_matrix() const {
			return new NCPA::BasicMatrix<T>();
		}
		virtual NCPA::Matrix<T>* build_matrix( size_t nrows, size_t ncols ) const {
			return new NCPA::BasicMatrix<T>( nrows, ncols );
		}

		// decompositions
		virtual NCPA::LUDecomposition<T>* build_lu_decomposition() const {
			return new NCPA::BasicLUDecomposition<T>();
		}
	};
}




#endif
