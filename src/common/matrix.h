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
#include "util.h"


namespace NCPA {

	// function templates for sorting multiple vectors by the same
	// criteria.  Stolen from
	// https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
	template <typename T, typename Compare>
	std::vector<std::size_t> sort_permutation(
	    const std::vector<T>& vec,
	    Compare compare) {
	    std::vector<std::size_t> p(vec.size());
	    std::iota(p.begin(), p.end(), 0);
	    std::sort(p.begin(), p.end(),
	        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
	    return p;
	}

	template <typename T>
	void apply_permutation_in_place(
	    std::vector<T>& vec,
	    const std::vector<std::size_t>& p) {
	    std::vector<bool> done(vec.size());
	    for (std::size_t i = 0; i < vec.size(); ++i) {
	        if (done[i]) {
	            continue;
	        }
	        done[i] = true;
	        std::size_t prev_j = i;
	        std::size_t j = p[i];
	        while (i != j) {
	            std::swap(vec[prev_j], vec[j]);
	            done[j] = true;
	            prev_j = j;
	            j = p[j];
	        }
	    }
	}

	template<typename T>
	class Matrix {

	public:
		virtual T get( size_t d1, size_t d2 ) const = 0;
		virtual T at( size_t d1, size_t d2 ) const = 0;
		virtual void set( size_t d1, size_t d2, T val ) = 0;
		virtual void ready() = 0;
		virtual bool is_ready() const = 0;
		virtual size_t rows() const = 0;
		virtual size_t columns() const = 0;
		virtual ~Matrix() { }

		virtual Matrix *add( const Matrix *second ) const = 0;
		virtual Matrix *multiply( const Matrix *second ) const = 0;
		virtual Matrix *transpose() const = 0;
		virtual void scale( T factor ) = 0;

		virtual void get_row( size_t row, size_t &n,
				size_t *&indices, T *&values ) const = 0;
		virtual void get_row( size_t row,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const = 0;
		virtual void get_row( size_t row, size_t &n,
				T *&values ) const = 0;
		virtual void get_row( size_t row,
				std::vector<T> &values ) const = 0;
		virtual void get_row_indices( size_t row, size_t &n,
				size_t *&indices ) const = 0;
		virtual void get_row_indices( size_t row,
				std::vector<size_t> &indices ) const = 0;
		virtual void get_column( size_t col, size_t &n,
				size_t *&indices, T *&values ) const = 0;
		virtual void get_column( size_t col,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const = 0;
		virtual void get_column( size_t col, size_t &n,
				T *&values ) const = 0;
		virtual void get_column( size_t col,
				std::vector<T> &values ) const = 0;
		virtual void get_column_indices( size_t col, size_t &n,
				size_t *&indices ) const = 0;
		virtual void get_column_indices( size_t col,
				std::vector<size_t> &indices ) const = 0;

		virtual void print() {
			for (size_t i = 0; i < rows(); i++) {
				std::cout << "[";
				for (size_t j = 0; j < columns(); j++) {
					std::cout << (j == 0 ? " " : ", ") << get(i,j);
				}
				std::cout << " ]" << std::endl;
			}
		}

	protected:
		virtual bool check_dimensions( size_t d1, size_t d2 ) const = 0;
	};




	template<typename T>
	class DenseMatrix : public Matrix<T> {

	protected:
		T **contents_;
		size_t rows_, cols_;
		bool ready_;

		virtual bool check_dimensions( size_t d1, size_t d2 ) const {
			return ( d1 < rows_ ) && ( d2 < cols_ );
		}

	public:
		DenseMatrix( size_t d1, size_t d2 ) {
			rows_ = (size_t)d1;
			cols_ = (size_t)d2;
			contents_ = NCPA::allocate_matrix<T>( rows_, cols_ );
			ready();
		}

		DenseMatrix( const NCPA::DenseMatrix<T> &in ) {
			rows_ = in.rows();
			cols_ = in.columns();
			contents_ = NCPA::allocate_matrix<T>( rows_, cols_ );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					contents_[i][j] = in.get( i, j );
				}
			}
			ready();
		}

		virtual ~DenseMatrix() {
			NCPA::free_matrix<T>( contents_, rows_, cols_ );
		}

		virtual T get( size_t d1, size_t d2 ) const {
			assert(ready_);
			check_dimensions( d1, d2 );
			return contents_[ d1 ][ d2 ];
		}

		virtual T at( size_t d1, size_t d2 ) const {
			return this->get( d1, d2 );
		}

		virtual void set( size_t d1, size_t d2, T val ) {
			check_dimensions( d1, d2 );
			contents_[ d1 ][ d2 ] = val;
		}

		// nothing special needs to be done to make ready
		virtual void ready() {
			ready_ = true;
		}

		virtual bool is_ready() const {
			return ready_;
		}

		virtual size_t rows() const {
			return rows_;
		}

		virtual size_t columns() const {
			return cols_;
		}

		virtual Matrix<T> *add( const Matrix<T> *term ) const {
			assert( this->rows() == term->rows() );
			assert( this->columns() == term->columns() );
			assert( ready_ );
			assert( term->is_ready() );

			Matrix<T> *sum = new DenseMatrix<T>( rows_, cols_ );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					sum->set( i, j, contents_[ i ][ j ] + term->get( i, j ) );
				}
			}

			sum->ready();
			return sum;
		}

		virtual Matrix<T> *multiply( const Matrix<T> *second ) const {
			assert( this->columns() == second->rows() );
			assert( ready_ );
			assert( second->is_ready() );

			Matrix<T> *product = new DenseMatrix<T>( rows_, second->columns() );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < second->columns(); j++) {
					T temp = 0;
					for (size_t k = 0; k < cols_; k++) {
						temp += contents_[i][k]
							* second->get(k,j);
					}
					product->set( i, j, temp );
				}
			}

			product->ready();
			return product;
		}

		virtual Matrix<T> *transpose() const {
			assert( ready_ );

			Matrix<T> *trans = new DenseMatrix<T>( cols_, rows_ );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					trans->set( j, i, contents_[ i ][ j ] );
				}
			}

			trans->ready();
			return trans;
		}

		virtual void scale( T factor ) {
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					contents_[ i ][ j ] *= factor;
				}
			}
		}

		virtual void get_row( size_t row, size_t &n,
				size_t *&indices, T *&values ) const {
			assert( row < rows_ );
			n = cols_;
			indices = NCPA::index_vector<size_t>( n );
			values  = NCPA::zeros<T>( n );
			std::memcpy( values, contents_[row], n*sizeof(T) );
		}

		virtual void get_row( size_t row, size_t &n,
				T *&values ) const {
			assert( row < rows_ );
			n = cols_;
			values  = NCPA::zeros<T>( n );
			std::memcpy( values, contents_[row], n*sizeof(T) );
		}

		virtual void get_row( size_t row,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const {
			assert( row < rows_ );
			indices.reserve( cols_ );
			values.reserve( cols_ );
			for (size_t i = 0; i < cols_; i++) {
				indices.push_back( i );
				values.push_back( contents_[row][i] );
			}
		}

		virtual void get_row( size_t row,
				std::vector<T> &values ) const {
			assert( row < rows_ );
			values.reserve( cols_ );
			for (size_t i = 0; i < cols_; i++) {
				values.push_back( contents_[row][i] );
			}
		}

		virtual void get_row_indices( size_t row, size_t &n,
				size_t *&indices ) const {
			assert( row < rows_ );
			n = cols_;
			indices = NCPA::index_vector<size_t>( n );
		}

		virtual void get_row_indices( size_t row,
				std::vector<size_t> &indices ) const {
			assert( row < rows_ );
			indices.reserve( cols_ );
			for (size_t i = 0; i < cols_; i++) {
				indices.push_back( i );
			}
		}


		virtual void get_column( size_t col, size_t &n,
				size_t *&indices, T *&values ) const {
			assert( col < cols_ );
			n = rows_;
			indices = NCPA::index_vector<size_t>( n );
			values  = NCPA::zeros<T>( n );
			for (size_t i = 0; i < n; i++) {
				values[ i ] = contents_[ i ][ col ];
			}
		}

		virtual void get_column( size_t col, size_t &n,
				T *&values ) const {
			assert( col < cols_ );
			n = rows_;
			values  = NCPA::zeros<T>( n );
			for (size_t i = 0; i < n; i++) {
				values[ i ] = contents_[ i ][ col ];
			}
		}

		virtual void get_column( size_t col,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const {
			assert( col < cols_ );
			indices.reserve( rows_ );
			values.reserve( rows_ );
			for (size_t i = 0; i < rows_; i++) {
				indices.push_back( i );
				values.push_back( contents_[i][col] );
			}
		}

		virtual void get_column( size_t col,
				std::vector<T> &values ) const {
			assert( col < cols_ );
			values.reserve( rows_ );
			for (size_t i = 0; i < rows_; i++) {
				values.push_back( contents_[i][col] );
			}
		}

		virtual void get_column_indices( size_t col, size_t &n,
				size_t *&indices ) const {
			assert( col < cols_ );
			n = rows_;
			indices = NCPA::index_vector<size_t>( n );
		}

		virtual void get_column_indices( size_t col,
				std::vector<size_t> &indices ) const {
			assert( col < cols_ );
			indices.reserve( rows_ );
			for (size_t i = 0; i < rows_; i++) {
				indices.push_back( i );
			}
		}

		virtual size_t get_row_capacity() const {
			return cols_;
		}


	};



	/****************************************************************
	 * SparseMatrix
	 ***************************************************************/

	template<typename T>
	class SparseMatrixRow {
	public:
		SparseMatrixRow() {}
		std::vector<size_t> *indices() {
			return &indices_;
		}
		std::vector<T> *values() {
			return &values_;
		}
		size_t size() const {
			assert(indices_.size() == values_.size());
			return values_.size();
		}
		void set( size_t index, T value ) {
			int ind = find_index( index );
			if (ind < 0) {
				indices_.push_back( index );
				values_.push_back( value );
			} else {
				values_[ ind ] = value;
			}
		}
		T get( size_t index ) {
			int ind = find_index( index );
			if (ind < 0) {
				return 0;
			} else {
				return values_[ ind ];
			}
		}
		void scale( T scalar ) {
			typename std::vector<T>::iterator it;
			for (it = values_.begin();
				it != values_.end(); ++it) {
				(*it) *= scalar;
			}
		}
		int find_index( size_t i ) {
			std::vector<size_t>::iterator it =
				std::find( indices_.begin(), indices_.end(), i );
			if (it == indices_.end()) {
				return -1;
			} else {
				return (int)(it - indices_.begin());
			}
		}
		void shrink_to_fit() {
			indices_.shrink_to_fit();
			values_.shrink_to_fit();
		}
		void sort() {
			std::vector<size_t> p = sort_permutation(indices_,
			    [](T const& a, T const& b){ return a < b; });
			apply_permutation_in_place( indices_, p );
			apply_permutation_in_place( values_, p );
		}
		void get_contents( size_t &n, size_t *&inds, T *&vals ) const {
			n = indices_.size();
			inds = new size_t[ n ];
			vals = new T[ n ];
			for (size_t i = 0; i < n; i++) {
				inds[ i ] = indices_[ i ];
				vals[ i ] = values_[ i ];
			}
		}

	protected:
		std::vector<size_t> indices_;
		std::vector<T>      values_;
	};

	template<typename T>
	class SparseMatrix : public Matrix<T> {

	public:
		SparseMatrix( size_t d1, size_t d2 ) {
			cols_ = d2;
			contents_.reserve( d1 );
			for (size_t i = 0; i < d1; i++) {
				SparseMatrixRow<T> *r = new SparseMatrixRow<T>();
				contents_.push_back( r );
			}
			ready_ = false;
		}

		SparseMatrix( const NCPA::Matrix<T> &in ) {
			cols_ = in.columns();
			contents_.reserve( in.rows() );
			size_t n;
			std::vector<size_t> indices;
			std::vector<T> values;
			for (size_t i = 0; i < in.rows(); i++) {
				SparseMatrixRow<T> *r = new SparseMatrixRow<T>();
				in.get_row( i, n, indices, values );
				for (size_t j = 0; j < n; j++) {
					r->set( indices[j], values[j] );
				}
			}
			ready();
		}


		~SparseMatrix() {
			typename std::vector< SparseMatrixRow<T> * >::iterator cit_;
			for (cit_ = contents_.begin(); cit_ != contents_.end(); ++cit_) {
				delete (*cit_);
			}
			contents_.clear();
		}

		virtual void set( size_t d1, size_t d2, T val ) {
			check_dimensions( d1, d2 );

			SparseMatrixRow<T> *row = contents_.at( d1 );
			row->set( d2, val );

			ready_ = false;
		}

		virtual T get( size_t d1, size_t d2 ) const {
			assert(ready_);
			check_dimensions( d1, d2 );
			return contents_.at( d1 )->get( d2 );
		}

		virtual T at( size_t d1, size_t d2 ) const {
			return this->get( d1, d2 );
		}

		virtual void ready() {
			typename std::vector< SparseMatrixRow<T> * >::iterator cit_;
			for (cit_ = contents_.begin(); cit_ != contents_.end(); ++cit_) {
				(*cit_)->shrink_to_fit();
				(*cit_)->sort();
			}
			ready_ = true;
		}

		virtual size_t rows() const {
			return contents_.size();
		}

		virtual size_t columns() const {
			return cols_;
		}

		virtual void get_row( size_t row, size_t &n,
				size_t *&indices, T *&values ) const {
			assert( row < rows() );
			assert( ready_ );
			contents_[ row ]->get_contents( n, indices, values );
		}

		virtual void get_row( size_t row, size_t &n,
				T *&values ) const {
			assert( row < rows() );
			assert( ready_ );
			size_t *indices;
			contents_[ row ]->get_contents( n, indices, values );
			delete [] indices;
		}

		virtual void get_row( size_t row,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const {
			assert( row < rows() );
			assert( ready_ );
			indices = *(contents_[ row ]->indices());
			values  = *(contents_[ row ]->values());
		}

		virtual void get_row( size_t row,
				std::vector<T> &values ) const {
			assert( row < rows() );
			assert( ready_ );
			values  = *(contents_[ row ]->values());
		}

		virtual void get_row_indices( size_t row, size_t &n,
				size_t *&indices ) const {
			assert( row < rows() );
			assert( ready_ );
			T *temp;
			contents_[row]->get_contents( n, indices, temp );
			delete [] temp;
		}

		virtual void get_row_indices( size_t row,
				std::vector<size_t> &indices ) const {
			assert( row < rows() );
			assert( ready_ );
			indices = *(contents_[ row ]->indices());
		}

		virtual void get_column( size_t col,
				std::vector<size_t> &inds,
				std::vector<T> &vals ) const {
			assert( col < columns() );
			assert( ready_ );
			inds.clear();
			vals.clear();
			for (size_t i = 0; i < contents_.size(); i++) {
				if (contents_[ i ]->find_index( col ) >= 0) {
					inds.push_back( i );
					vals.push_back( contents_[ i ]->get( i ) );
				}
			}
		}

		virtual void get_column( size_t col,
				std::vector<T> &vals ) const {
			assert( col < columns() );
			assert( ready_ );
			vals.clear();
			for (size_t i = 0; i < contents_.size(); i++) {
				if (contents_[ i ]->find_index( col ) >= 0) {
					vals.push_back( contents_[ i ]->get( i ) );
				}
			}
		}

		virtual void get_column( size_t col, size_t &n,
				size_t *&indices, T *&values ) const {
			assert( col < cols_ );
			assert( ready_ );
			std::vector< size_t > inds;
			std::vector< T > vals;
			get_column( col, inds, vals );
			n = inds.size();
			indices = new size_t[ n ];
			values  = new T[ n ];
			for (size_t i = 0; i < n; i++) {
				indices[ i ] = inds[ i ];
				values[ i ] = vals[ i ];
			}
		}

		virtual void get_column( size_t col, size_t &n,
				T *&values ) const {
			assert( col < cols_ );
			assert( ready_ );
			std::vector< size_t > inds;
			std::vector< T > vals;
			get_column( col, inds, vals );
			n = inds.size();
			values  = new T[ n ];
			for (size_t i = 0; i < n; i++) {
				values[ i ] = vals[ i ];
			}
		}

		virtual void get_column_indices( size_t col,
				std::vector<size_t> &inds ) const {
			assert( col < columns() );
			assert( ready_ );
			inds.clear();
			for (size_t i = 0; i < contents_.size(); i++) {
				if (contents_[ i ]->find_index( col ) >= 0) {
					inds.push_back( i );
				}
			}
		}

		virtual void get_column_indices( size_t col, size_t &n,
				size_t *&indices ) const {
			assert( col < cols_ );
			assert( ready_ );
			std::vector< size_t > inds;
			get_column( col, inds );
			n = inds.size();
			indices = new size_t[ n ];
			for (size_t i = 0; i < n; i++) {
				indices[ i ] = inds[ i ];
			}
		}

		virtual Matrix<T> *add( const Matrix<T> *term ) const {
			assert( rows() == term->rows() );
			assert( columns() == term->columns() );
			assert( is_ready() );
			assert( term->is_ready() );

			Matrix<T> *sum = new SparseMatrix<T>( rows(), columns() );
			for (size_t i = 0; i < rows(); i++) {
				std::vector<size_t> u( rows() + term->rows() ), inds1, inds2;
				std::vector<size_t>::iterator it;
				this->get_row( i, inds1 );
				term->get_row( i, inds2 );
				it = std::set_union(
					inds1.begin(), inds1.end(),
					inds2.begin(), inds2.end(),
					u.begin() );
				size_t n = it - u.begin();
				for (it = u.begin(); it != u.end(); ++it) {
					sum->set( i, *it,
						this->get( i, *it ) + term->get( i, *it ) );
				}
			}

			sum->ready();
			return sum;
		}

		virtual Matrix<T> *multiply( const Matrix<T> *term ) const {
			assert( columns() == term->rows() );
			assert( is_ready() );
			assert( term->is_ready() );

			Matrix<T> *product = new SparseMatrix<T>(
				rows(), term->columns() );

			for (size_t i = 0; i < rows(); i++) {
				for (size_t j = 0; j < term->columns(); j++) {
					std::vector<size_t> u( rows() + term->rows() ),
							inds1, inds2;
					std::vector<size_t>::iterator it;
					this->get_row( i, inds1 );
					term->get_column( j, inds2 );
					it = std::set_intersection(
						inds1.begin(), inds1.end(),
						inds2.begin(), inds2.end(),
						u.begin() );
					size_t n = it - u.begin();
					T factor = 0;
					// std::cout << "product[" << i << "][" << j << "] = ";
					for (size_t k = 0; k < n; k++) {
						// if (k > 0) {
						// 	std::cout << " + ";
						// }
						T f1 = this->get( i, u[k] );
						T f2 = term->get( u[k], j );
						factor += f1*f2;
						// std::cout << f1 << "*" << f2;

					}
					// std::cout << std::endl;
					product->set( i, j, factor );
				}
			}

			product->ready();
			return product;
		}

		virtual Matrix<T> *transpose() const {
			assert( ready_ );

			Matrix<T> *trans = new SparseMatrix<T>( columns(), rows() );
			for (size_t i = 0; i < rows(); i++) {
				std::vector<size_t> inds = *(contents_[ i ]->indices());
				for (std::vector<size_t>::iterator j = inds.begin();
						j != inds.end(); ++j) {
					trans->set( *j, i, get( i, *j ) );
				}
			}

			trans->ready();
			return trans;
		}

		virtual void scale( T scalar ) {
			typename std::vector< SparseMatrixRow<T> * >::iterator cit_;
			for (cit_ = contents_.begin(); cit_ != contents_.end(); ++cit_) {
				(*cit_)->scale( scalar );
			}
		}

		virtual bool is_ready() const {
			return ready_;
		}


	protected:
		std::vector< SparseMatrixRow<T> * > contents_;
		// std::vector< SparseMatrixRow<T> * >::iterator cit_;
		size_t cols_;
		bool ready_;

		virtual bool check_dimensions( size_t d1, size_t d2 ) const {
			assert( d1 < rows() );
			assert( d2 < columns() );
		}
	};

	template<typename T>
	class DiagonalMatrix : public SparseMatrix<T> {

	public:
		DiagonalMatrix( size_t d1, size_t d2, T val )
				: SparseMatrix<T>( d1, d2, 1 ) {
			for (size_t i = 0; i < NCPA::min(d1,d2); i++) {
				set( i, i, val );
			}
			this->ready_ = true;
			this->contiguous_ = true;
		}

		DiagonalMatrix( size_t d1, size_t d2, T *vals )
				: SparseMatrix<T>( d1, d2, 1 ) {
			for (size_t i = 0; i < NCPA::min(d1,d2); i++) {
				set( i, i, vals[i] );
			}
			this->ready_ = true;
			this->contiguous_ = true;
		}
	};

	template<typename T>
	class TridiagonalMatrix : public SparseMatrix<T> {

	public:
		TridiagonalMatrix( size_t d1, size_t d2, T diagval,
			T offdiagval_lower, T offdiagval_upper )
			: SparseMatrix<T>( d1, d2, 3 ) {

			size_t minind = NCPA::min(d1,d2);
			for (size_t i = 0; i < minind; i++) {
				if (i > 0) {
					set( i, i-1, offdiagval_lower );
				}
				set( i, i, diagval );
				if (i < (minind-1)) {
					set( i, i+1, offdiagval_upper );
				}
			}
			this->ready_ = true;
			this->contiguous_ = true;
		}

		TridiagonalMatrix( size_t d1, size_t d2, T *diagvals,
			T *offdiagvals_lower, T *offdiagvals_upper )
			: SparseMatrix<T>( d1, d2, 3 ) {

			size_t minind = NCPA::min(d1,d2);
			for (size_t i = 0; i < minind; i++) {
				if (i > 0) {
					set( i, i-1, offdiagvals_lower[i] );
				}
				set( i, i, diagvals[i] );
				if (i < (minind-1)) {
					set( i, i+1, offdiagvals_upper[i] );
				}
			}
			this->ready_ = true;
			this->contiguous_ = true;

		}
	};



}







#endif