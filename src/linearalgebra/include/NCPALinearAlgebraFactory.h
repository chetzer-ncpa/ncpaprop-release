#ifndef NCPA_LINEARALGEBRA__LINEARALGEBRAFACTORY_H_INCLUDED
#define NCPA_LINEARALGEBRA__LINEARALGEBRAFACTORY_H_INCLUDED

#include "NCPAVector.h"
#include "NCPAMatrix.h"
#include "NCPALUDecomposition.h"

namespace NCPA {
	template<typename T>
	class LinearAlgebraFactory  {
	public:
		virtual ~LinearAlgebraFactory() {}

		// vectors
		virtual NCPA::Vector<T>* build_vector() const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t nvals ) const = 0;
		virtual NCPA::Vector<T>* build_vector( std::vector<T> vals ) const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::vector<size_t> inds, std::vector<T> vals ) const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t n,
				size_t nvals, const size_t *inds, const T *vals ) const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::initializer_list<size_t> inds,
				std::initializer_list<T> vals ) const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::vector<size_t> inds, T val ) const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t n, size_t nvals,
				const size_t *inds, T val ) const = 0;
		virtual NCPA::Vector<T>* build_vector( size_t n,
				std::initializer_list<size_t> inds, T val ) const = 0;

		// matrices
		virtual NCPA::Matrix<T>* build_matrix() const = 0;
		virtual NCPA::Matrix<T>* build_matrix( size_t nrows, size_t ncols ) const = 0;

		// decompositions
		virtual NCPA::LUDecomposition<T>* build_lu_decomposition() const = 0;
	};
}

#endif
