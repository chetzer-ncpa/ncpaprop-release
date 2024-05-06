#ifndef NCPA__LINEARALGEBRA_SPARSEVECTOR_H_INCLUDED_
#define NCPA__LINEARALGEBRA_SPARSEVECTOR_H_INCLUDED_

#include <unordered_map>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace NCPA { template<typename T> class SparseVector; }

template<typename T>
void swap( NCPA::SparseVector<T> &a, NCPA::SparseVector<T> &b ) noexcept;

namespace NCPA {

	template<typename T>
	class SparseVector : public std::unordered_map<size_t,T> {

	public:
		SparseVector() = default;

		SparseVector( size_t n ) : std::unordered_map<size_t,T>() {
			this->reserve(n);
		}

		SparseVector( std::vector<size_t> inds, std::vector<T> vals ) :
			NCPA::SparseVector<T>( inds.size() ) {
			if (inds.size() != vals.size()) {
				throw std::out_of_range( "Index and value vectors are of different sizes!" );
			}
			for (size_t i = 0; i < inds.size(); i++) {
				(*this)[inds[i]] = vals[i];
			}
		}

		SparseVector( size_t n, const size_t *inds, const T *vals ) :
			NCPA::SparseVector<T>(
					std::vector<size_t>( inds, inds+n ),
					std::vector<T>( vals, vals+n ) )
			{}

		SparseVector( std::initializer_list<size_t> inds, std::initializer_list<T> vals) :
			NCPA::SparseVector<T>( std::vector<size_t>( inds ), std::vector<T>( vals ) )
			{}

		SparseVector( std::vector<size_t> inds, T val ) :
			NCPA::SparseVector<T>( inds, std::vector<T>( inds.size(), val ) )
			{}

		SparseVector( size_t n, const size_t *inds, T val ) :
			NCPA::SparseVector<T>( std::vector<size_t>( inds, inds+n ),
				std::vector<T>( n, val ) )
			{}

		SparseVector( std::initializer_list<size_t> inds, T val ) :
			NCPA::SparseVector<T>( std::vector<size_t>( inds ),
					std::vector<T>( inds.size(), val ) )
			{}
	};


}

#endif
