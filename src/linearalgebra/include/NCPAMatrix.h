#ifndef NCPA_MATRIX_H_INCLUDED
#define NCPA_MATRIX_H_INCLUDED

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "NCPACommon.h"
#include "NCPASparseVector.h"

namespace NCPA { template<typename T> class Matrix; }

template<typename T>
void swap( NCPA::Matrix<T> &a, NCPA::Matrix<T> &b ) noexcept {}

template<typename T>
bool operator==(const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b);
template<typename T>
bool operator!=(const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b);

namespace NCPA {

	template<typename T>
	class Matrix {

	public:
		Matrix() {}
		Matrix( const Matrix &other ) {}
		Matrix( Matrix &&other ) noexcept {}
		virtual ~Matrix() {}
		friend void swap<T>( Matrix<T> &a, Matrix<T> &b ) noexcept;

		// Interface
		virtual Matrix<T>* clone() const = 0;
		virtual NCPA::Matrix<T> * initialize( size_t nrows, size_t ncols ) = 0;

		// element access
		virtual size_t rows() const = 0;
		virtual size_t columns() const = 0;
		virtual T get( size_t d1, size_t d2 ) const = 0;
		virtual Matrix<T>* set( size_t d1, size_t d2, T val ) = 0;


		// row and column access
		virtual void get_row( size_t rownum, size_t &n, T* &values ) const = 0;
		virtual void get_row( size_t rownum, std::vector<T> &values ) const = 0;
		virtual void get_row( size_t rownum, NCPA::SparseVector<T> &row )  const = 0;
		virtual void get_column( size_t colnum, size_t &n, T* &values ) const = 0;
		virtual void get_column( size_t colnum, std::vector<T> &values ) const = 0;
		virtual void get_column( size_t colnum, NCPA::SparseVector<T> &row ) const = 0;
//		virtual std::vector<T> get_diagonal() const = 0;
//		virtual void get_diagonal( size_t &n, T* &vals ) const = 0;
		virtual std::vector<size_t> get_column_indices(size_t row) const = 0;

		virtual Matrix<T>* set_row( size_t rownum, size_t nvals,
				const size_t* column_indices, const T* values ) = 0;
		virtual Matrix<T>* set_row( size_t rownum, const std::vector<T> &values ) = 0;
		virtual Matrix<T>* set_row( size_t rownum,
				const std::vector<size_t> &column_indices,
				const std::vector<T> &values ) = 0;
		virtual Matrix<T>* set_row( size_t rownum,
				const NCPA::SparseVector<T> &values ) = 0;
		virtual Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) = 0;
		virtual Matrix<T>* set_column( size_t rownum, const std::vector<T> &values ) = 0;
		virtual Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> &row_indices,
				const std::vector<T> &values ) = 0;
		virtual Matrix<T>* set_column( size_t colnum,
				const NCPA::SparseVector<T> &values ) = 0;

		// status
		virtual Matrix<T>* ready() = 0;
		virtual bool is_ready() const = 0;

		// operations
		virtual Matrix<T> *add( const NCPA::Matrix<T> *second ) = 0;
		virtual Matrix<T> *subtract( const NCPA::Matrix<T> *second ) = 0;
		virtual Matrix<T> *transpose() = 0;
		virtual Matrix<T> *scale( T factor ) = 0;
		virtual Matrix<T> *negative() {
			this->scale(-1.0);
			return this;
		}


		// concrete base methods below here.  Base classes may implement their own
		// versions
		static void multiply(
				const NCPA::Matrix<T> *first,
				const NCPA::Matrix<T> *second,
				NCPA::Matrix<T> *&product ) {
			if (first->columns() != second->rows()) {
				std::ostringstream oss;
				oss << "Size mismatch in matrix product: Columns in first matrix ("
						<< first->columns() << ") must equal rows in second matrix ("
						<< second->rows() << ")";
				throw std::out_of_range( oss.str() );
			}
			if (product == nullptr) {
				throw std::logic_error( "No product matrix provided" );
			}
			product->initialize(first->rows(), second->columns());
			for (size_t r = 0; r < first->rows(); r++) {
				for (size_t c = 0; c < second->columns(); c++) {
					T prod = 0;
					for (size_t i = 0; i < first->columns(); i++) {
						prod += first->get(r,i) * second->get(i,c);
					}
					product->set( r, c, prod );
				}
			}
		}

		virtual std::vector<T> get_diagonal() const {
			size_t n = rows() < columns() ? rows() : columns();
			std::vector<T> diag(n);
			for (size_t ii = 0; ii < n; ii++) {
				diag[ii] = this->get(ii,ii);
			}
			return diag;
		}

		virtual void get_diagonal( size_t &n, T* &diag ) const {
			size_t nn = rows() < columns() ? rows() : columns();
			if (diag == nullptr) {
				n = nn;
				diag = NCPA::zeros<T>( n );
			} else {
				if (n == 0) {
					n = nn;
				} else if (n != nn) {
					std::ostringstream oss;
					oss << "Matrix diagonal size mismatch: diagonal has "
							<< nn << " elements, but " << n << " requested.";
					throw std::out_of_range( oss.str() );
				}
			}
			for (size_t ii = 0; ii < n; ii++) {
				diag[ii] = this->get(ii,ii);
			}
		}

		virtual void get_offdiagonal( int offset, size_t &n, T* &values ) const {
			std::vector<T> v = this->get_offdiagonal( offset );
			if (values == nullptr) {
				n = v.size();
				values = NCPA::zeros<T>( n );
			} else {
				if (n == 0) {
					n = v.size();
				} else if (n != v.size()) {
					std::ostringstream oss;
					oss << "Matrix offdiagonal size mismatch: offdiagonal has "
							<< v.size() << " elements, but " << n << " requested.";
					throw std::out_of_range( oss.str() );
				}
			}
			std::copy( v.cbegin(), v.cend(), values );
		}

		virtual std::vector<T> get_offdiagonal( int offset ) const {
			if (offset == 0) {
				return this->get_diagonal();
			}
			int ndiag;
			size_t rstart, cstart;
			this->calculate_offdiagonal_parameters_( offset, ndiag, rstart, cstart );
			std::vector<T> values(ndiag);
			for (size_t i = 0; i < (size_t)ndiag; i++) {
				values[i] = this->get( rstart+i, cstart+i );
			}
			return values;
		}

		virtual NCPA::Matrix<T>* set_diagonal( size_t n, const T* values ) {
			this->check_dimensions( n, n );
			for (size_t i = 0; i < n; i++) {
				this->set( i, i, values[i] );
			}
			return this;
		}

		virtual NCPA::Matrix<T>* set_diagonal( const std::vector<T> &values ) {
			this->check_dimensions( values.size(), values.size() );
			for (size_t i = 0; i < values.size(); i++) {
				this->set( i, i, values[i] );
			}
			return this;
		}

		virtual NCPA::Matrix<T>* set_diagonal( T val ) {
			size_t n = this->calculate_diagonal_length_();
			for (size_t i = 0; i < n; i++) {
				this->set( i, i, val );
			}
			return this;
		}

		virtual NCPA::Matrix<T>* set_offdiagonal( int offset, T val ) {
			int ndiag;
			size_t rstart, cstart;
			calculate_offdiagonal_parameters_( offset, ndiag, rstart, cstart );
			for (size_t i = 0; i < ndiag; i++) {
				this->set( rstart+i, cstart+i, val );
			}
			return this;
		}

		virtual NCPA::Matrix<T>* set_offdiagonal( int offset,
				size_t n, const T* values ) {
			int ndiag;
			size_t rstart, cstart;
			calculate_offdiagonal_parameters_( offset, ndiag, rstart, cstart );
			if (n != ndiag) {
				std::ostringstream oss;
				oss << "Size mismatch: off-diagonal " << offset << " has "
						<< ndiag << " elements for matrix size " << rows()
						<< "x" << columns() << ", but " << n << " elements provided.";
				throw std::out_of_range(oss.str());
			}
			for (size_t i = 0; i < n; i++) {
				this->set( rstart+i, cstart+i, values[i] );
			}
			return this;
		}

		virtual NCPA::Matrix<T>* set_offdiagonal( int offset,
				const std::vector<T> &values ) {
			int ndiag;
			size_t rstart, cstart;
			calculate_offdiagonal_parameters_( offset, ndiag, rstart, cstart );
			if (values.size() != ndiag) {
				std::ostringstream oss;
				oss << "Size mismatch: off-diagonal " << offset << " has "
						<< ndiag << " elements for matrix size " << rows()
						<< "x" << columns() << ", but " << values.size()
						<< " elements provided.";
				throw std::out_of_range(oss.str());
			}
			for (size_t i = 0; i < values.size(); i++) {
				this->set( rstart+i, cstart+i, values[i] );
			}
			return this;
		}


		virtual Matrix<T>* set_row( size_t rownum,
						std::initializer_list<size_t> column_indices,
						std::initializer_list<T> values ) {
			return this->set_row( rownum,
					std::vector<size_t>(column_indices),
					std::vector<T>( values ) );
		}
		virtual Matrix<T>* set_column( size_t colnum,
				std::initializer_list<size_t> row_indices,
				std::initializer_list<T> values ) {
			return this->set_column( colnum,
					std::vector<size_t>(row_indices),
					std::vector<T>( values ) );
		}
//
//		NCPA::SparseVector<T>& operator[]( int rownum ) {
//			NCPA::SparseVector<T> row;
//			this->get_row( rownum, row );
//			return row;
//		}

		// operators
		friend bool operator!=( const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b ) {
			return a.is_unequal_( b );
		}
		friend bool operator==( const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b ) {
			return !(a.is_unequal_( b ));
		}


		virtual void check_dimensions( size_t d1, size_t d2 ) const {
			this->check_rows_(d1);
			this->check_columns_(d2);
		}

		virtual void print( std::ostream &o = std::cout ) {
			o << "[" << std::endl;
			for (size_t r = 0; r < rows(); r++) {
				o << "  [ ";
				for (size_t c = 0; c < columns(); c++) {
					o << get(r,c);
					if (c != (columns()-1)) {
						o << ", ";
					} else {
						o << " ]" << std::endl;
					}
				}
			}
			o << "]" << std::endl;
		}

	protected:
		virtual bool is_unequal_( const Matrix<T> &other ) const = 0;

		virtual void check_rows_( size_t d1 ) const {
			if (d1 > this->rows()) {
				std::ostringstream oss;
				oss << "Requested row " << d1
						<< " is out of range for Matrix with "
						<< rows() << " rows.";
				throw std::out_of_range( oss.str() );
			}
		}
		virtual void check_columns_( size_t d1 ) const {
			if (d1 > this->columns()) {
				std::ostringstream oss;
				oss << "Requested column " << d1
						<< " is out of range for Matrix with "
						<< columns() << " columns.";
				throw std::out_of_range( oss.str() );
			}
		}
		virtual void check_sizes_match_( const Matrix<T> *other ) const {
			if (this->rows() != other->rows() || this->columns() != other->columns()) {
				std::ostringstream oss;
				oss << "Matrix size mismatch between RHS matrix of size "
						<< other->rows() << "x" << other->columns()
						<< " and LHS matrix of size " << this->rows() << "x"
						<< this->columns() << std::endl;
				throw std::out_of_range( oss.str() );
			}
		}
		virtual void check_sizes_match_transpose_( const Matrix<T> *other ) const {
			if (this->rows() != other->columns() || this->columns() != other->rows()) {
				std::ostringstream oss;
				oss << "Matrix transpose size mismatch between RHS matrix of size "
						<< other->columns() << "x" << other->rows()
						<< " and LHS matrix of size " << this->rows() << "x"
						<< this->columns() << std::endl;
				throw std::out_of_range( oss.str() );
			}
		}
		virtual size_t calculate_diagonal_length_() const {
			return NCPA::min<size_t>( this->rows(), this->columns() );
		}
		virtual void calculate_offdiagonal_parameters_( int offset,
				int &ndiag, size_t &rstart, size_t &cstart ) const {
			int ndiag_max = (int)(this->calculate_diagonal_length_());
			int sizediff = this->columns() - this->rows();
			rstart = 0;
			cstart = 0;

			// square matrix: ndiag = rows() - abs(offset)
			if (sizediff == 0) {
				ndiag = (int)(this->rows()) - std::abs(offset);
			} else if (NCPA::sign<int>(offset) == NCPA::sign<int>(sizediff)) {
				if (std::abs(offset) <= std::abs(sizediff)) {
					ndiag = ndiag_max;
				} else {
					ndiag = ndiag_max + std::abs(sizediff) - std::abs(offset);
				}
			} else {
				ndiag = ndiag_max - std::abs(offset);
			}
			if (ndiag <= 0) {
				std::ostringstream oss;
				oss << "Offdiagonal " << offset << " is too large for a matrix of size "
						<< this->rows() << "x" << this->columns() << ".";
				throw std::out_of_range(oss.str());
			}

			if (offset < 0) {
				rstart = (size_t)std::abs(offset);
			} else {
				cstart = (size_t)std::abs(offset);
			}
		}
	};
}
//
//namespace NCPA {
//
//	// function templates for sorting multiple vectors by the same
//	// criteria.  Stolen from
//	// https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
//	template <typename T, typename Compare>
//	std::vector<std::size_t> sort_permutation(
//	    const std::vector<T>& vec,
//	    Compare compare) {
//	    std::vector<std::size_t> p(vec.size());
//	    std::iota(p.begin(), p.end(), 0);
//	    std::sort(p.begin(), p.end(),
//	        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
//	    return p;
//	}
//
//	template <typename T>
//	void apply_permutation_in_place(
//	    std::vector<T>& vec,
//	    const std::vector<std::size_t>& p) {
//	    std::vector<bool> done(vec.size());
//	    for (std::size_t i = 0; i < vec.size(); ++i) {
//	        if (done[i]) {
//	            continue;
//	        }
//	        done[i] = true;
//	        std::size_t prev_j = i;
//	        std::size_t j = p[i];
//	        while (i != j) {
//	            std::swap(vec[prev_j], vec[j]);
//	            done[j] = true;
//	            prev_j = j;
//	            j = p[j];
//	        }
//	    }
//	}
//}


#endif
