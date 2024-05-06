/*
NCPABasicMatrix.h

Provides a concrete instance of a sparse matrix, where most values are expected to be
zero and so only the nonzero elements are stored.  The extra overhead involved in
identifying the nonzero elements is compensated for by many fewer operations in
matrix arithmetic and solution operations.  Specific types of sparse matrices, such
as diagonal and tridiagonal matrices, have their own classes derived from this class.

Provided under the following license:

All rights reserved.
Copyright 2016 University of Mississippi

Developed by: National Center for Physical Acoustics, developer: Roger Waxler

Permission is hereby granted by the University of Mississippi (Grantor), free of charge, to any person (Grantee) obtaining a
copy of this software and associated documentation files (the ”Software”), to deal with the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
− Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
− Any modifications to the source code must carry prominent notices stating that the files were changed.
− Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following
disclaimers in the documentation and/or other materials provided with the distribution;
− Neither the names of the University of Mississippi, the National Center for Physical Acoustics, nor the names of its
developers may be used to endorse or promote products derived from this Software without specific prior written permission;
− Grantee must give any other recipients of the Work or Derivative Works a copy of this License;

THE SOFTWARE IS PROVIDED ”AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS WITH THE SOFTWARE.
*/

/**
 * @file NCPABasicMatrix.h
 * @author Claus Hetzer
 * @version 1.0.0
 * @date 2023-11-10
 * @brief Provides a template class for sparse matrices.
 */

#ifndef NCPA__LINEARALGEBRA_NCPABASICMATRIX_H_INCLUDED_
#define NCPA__LINEARALGEBRA_NCPABASICMATRIX_H_INCLUDED_

#include "NCPACommon.h"
#include <map>
#include <stdexcept>
#include <sstream>
#include <vector>
#include "NCPAMatrix.h"
#include "NCPABasicVector.h"

namespace NCPA { template<typename T> class BasicMatrix; }

template<typename T>
void swap( NCPA::BasicMatrix<T> &a, NCPA::BasicMatrix<T> &b ) noexcept;
template<typename T>
bool operator==(const NCPA::BasicMatrix<T> &a, const NCPA::BasicMatrix<T> &b);
template<typename T>
bool operator!=(const NCPA::BasicMatrix<T> &a, const NCPA::BasicMatrix<T> &b);

namespace NCPA {
	template<typename T>
	class BasicMatrix : public NCPA::Matrix<T> {
	public:


		// constructors
		BasicMatrix() : NCPA::Matrix<T>() {}
		BasicMatrix( size_t nrows, size_t ncols ) : NCPA::Matrix<T>() {
			this->initialize( nrows, ncols );
		}
		BasicMatrix( const NCPA::BasicMatrix<T> &other ) : NCPA::Matrix<T>(other) {
			for (size_t i = 0; i < other.rows_.size(); i++) {
				this->rows_.push_back( new NCPA::BasicVector<T>( *(other.rows_[i]) ) );
			}
			for (auto it = other.cols_cache_.cbegin();
					it != other.cols_cache_.cend(); ++it) {
				this->cols_cache_[it->first] =
						new NCPA::BasicVector<T>( *(it->second) );
			}
//			this->rows_ = other.rows_;
			this->ncols_ = other.ncols_;

		}
		BasicMatrix( NCPA::BasicMatrix<T> &&other ) noexcept : NCPA::Matrix<T>() {
			::swap(*this,other);
		}
		virtual ~BasicMatrix() {
			this->clear();
		}

		//swap
		friend void ::swap<T>( NCPA::BasicMatrix<T> &a, NCPA::BasicMatrix<T> &b ) noexcept;

		// assignment operator
		NCPA::BasicMatrix<T>& operator=( NCPA::BasicMatrix<T> other ) {
			swap( *this, other );
			return *this;
		}

		virtual NCPA::Matrix<T> * clone() const {
			return static_cast<NCPA::Matrix<T> *>(
					new NCPA::BasicMatrix<T>( *this ) );
		}

		virtual size_t rows() const {
			return rows_.size();
		}

		virtual size_t columns() const {
			return ncols_;
		}

		virtual NCPA::Matrix<T> * clear() {
			for (auto it = this->rows_.begin(); it != this->rows_.end(); ++it) {
				if (*it != nullptr) {
					delete *it;
					*it = nullptr;
				}
			}
			this->rows_.clear();
			this->clear_cache();
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		// initializer
		virtual NCPA::Matrix<T> * initialize( size_t nrows, size_t ncols ) {
			this->assert_not_finalized();
			this->clear();
			rows_ = std::vector<NCPA::BasicVector<T>*>( nrows, nullptr );
			ncols_ = ncols;
			for (auto it = rows_.begin(); it != rows_.end(); ++it) {
				*it = new NCPA::BasicVector<T>( ncols );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> * finalize() {
			// @todo Populate column cache?
			this->finalized_ = true;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> * clear_cache() {
			for (auto it = this->cols_cache_.begin();
					it != this->cols_cache_.end(); ++it) {
				if (it->second != nullptr) {
					delete it->second;
				}
			}
			this->cols_cache_.clear();
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual T get( size_t row, size_t col ) const {
			this->assert_element_in_range(row,col);
			return rows_[row]->get(col);
//			auto it = rows_[row].find(col);
//			if (it == rows_[row].cend()) {
//				return 0;
//			} else {
//				return it->second;
//			}
		}

		virtual NCPA::Matrix<T>* set( size_t row, size_t col, T val ) {
			this->assert_not_finalized();
			this->assert_element_in_range( row, col );
			rows_[row]->set( col, val );
//			rows_[row][col] = val;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		/**
		 * Zero out a row of the matrix.
		 * @param[in] rownum The number of the row to clear, indexed to 0.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* clear_row( size_t rownum ) {
			this->assert_row_in_range( rownum );
			rows_[rownum]->clear()->resize( ncols_ );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual void get_row( size_t rownum, size_t &n, T* &values ) const {
			this->assert_row_in_range( rownum );
			if (values == nullptr) {
				n = columns();
				values = NCPA::zeros<T>( n );
			} else if (n != columns()) {
				std::ostringstream oss;
				oss << "Matrix row size mismatch: Matrix has " << rows()
						<< " rows, but " << n << " values requested.";
				throw std::out_of_range( oss.str() );
			}
			rows_[rownum]->as_array( n, values );
//			for (auto it = rows_[rownum].cbegin(); it != rows_[rownum].cend(); ++it) {
//				values[ it->first ] = it->second;
//			}
		}

		virtual void get_row( size_t rownum, std::vector<T> &values ) const {
			this->assert_row_in_range( rownum );
			values = rows_[rownum]->as_std();
//			values.resize( columns() );
//
//			for (auto it = rows_[rownum].cbegin(); it != rows_[rownum].cend(); ++it) {
//				values[ it->first ] = it->second;
//			}
		}

		virtual const NCPA::Vector<T>* get_row( size_t rownum ) const {
			this->assert_row_in_range( rownum );
			return rows_[rownum];
//			NCPA::Vector<T> *v = rows_[rownum].clone();
//			this->get_row( rownum, v );
//			return v;
		}

		virtual void get_row( size_t rownum, NCPA::Vector<T> *values ) const {
			if (values == nullptr) {
				throw std::invalid_argument( "Vector pointer is null!");
			}
			this->assert_row_in_range( rownum );
			values->copy( rows_[rownum] );
		}

		/**
		 * Zero out a column of the matrix.
		 * @param[in] colnum The number of the column to clear, indexed to 0.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* clear_column( size_t colnum ) {
			this->assert_column_in_range( colnum );
			for (auto it = rows_.begin(); it != rows_.end(); ++it) {
				(*it)->clear( colnum );
			}
			auto cache_it = cols_cache_.find( colnum );
			if (cache_it != cols_cache_.end()) {
				cols_cache_.erase( cache_it );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual const NCPA::Vector<T>* get_column( size_t colnum ) const {
			this->assert_column_in_range( colnum );
			auto it = this->cols_cache_.find(colnum);
			if (it != this->cols_cache_.cend()) {
				return it->second;
			}

			NCPA::BasicVector<T> *v = new NCPA::BasicVector<T>( this->rows() );
			for (size_t i = 0; i < this->rows(); i++) {
				std::vector<size_t> inds = rows_[i]->get_defined_indices();
				std::vector<size_t>::iterator it = std::find(
						inds.begin(), inds.end(), colnum );
				if (it != inds.end()) {
					v->set( i, rows_[i]->get(colnum) );
				}
			}
			cols_cache_[ colnum ] = v;
			return v;
//
////			NCPA::Vector<T> *v = NCPA::Vector<T>::build(NCPA::vector_t::SPARSE);
//			NCPA::Vector<T> *v = nullptr;
//			this->get_column( colnum, v );
//			this->cols_cache_[ colnum ] = v;
//			return v;
		}

		virtual void get_column( size_t colnum, size_t &n, T* &values ) const {
			this->assert_column_in_range( colnum );
			if (values == nullptr) {
				n = rows();
				values = NCPA::zeros<T>( n );
			} else if (n != rows()) {
				std::ostringstream oss;
				oss << "Matrix column size mismatch: Matrix has " << columns()
						<< " columns, but " << n << " values requested.";
				throw std::out_of_range( oss.str() );
			}
			this->get_column( colnum )->as_array( n, values );
//
//
//			this->assert_column_in_range( colnum );
//			auto it = this->cols_cache_.find(colnum);
//			if (it != this->cols_cache_.cend()) {
//				it->as_array( n, values );
//			}
//
//
//			size_t ii = 0;
//			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
//				if (it->find(colnum) != it->cend()) {
//					values[ii] = it->at(colnum);
//				} else {
//					values[ii] = 0.0;
//				}
//				ii++;
//			}
		}

		virtual void get_column( size_t colnum, std::vector<T> &values ) const {
			this->assert_column_in_range( colnum );
			values = this->get_column( colnum )->as_std();
//			values.resize( rows() );
//			size_t ii = 0;
//			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
//				if (it->find(colnum) != it->cend()) {
//					values[ii] = it->at(colnum);
//				} else {
//					values[ii] = 0.0;
//				}
//				ii++;
//			}
		}

		virtual void get_column( size_t colnum, NCPA::Vector<T> *values ) const {
			if (values == nullptr) {
				throw std::invalid_argument( "Vector pointer is null!");
			}
			this->assert_column_in_range( colnum );
			values->copy( this->get_column( colnum ) );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum, size_t nvals,
						const size_t * column_indices, const T* values ) {
			this->assert_not_finalized();
			this->assert_row_in_range( rownum );
			rows_[rownum]->clear()->resize(this->ncols_);
			for (size_t i = 0; i < nvals; i++) {
				this->assert_column_in_range( column_indices[i] );
				rows_[rownum]->set(column_indices[i], values[i]);
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<T> &values ) {
			this->assert_not_finalized();
			this->assert_row_in_range( rownum );
			this->assert_columns_match( values.size() );
			rows_[rownum]->from_std( values );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<size_t> &column_indices,
						const std::vector<T> &values ) {
			this->assert_not_finalized();
			this->assert_row_in_range( rownum );
			rows_[rownum]->clear()->resize(this->ncols_);
			for (size_t i = 0; i < values.size(); i++) {
				this->assert_column_in_range( column_indices[i] );
				rows_[rownum]->set(column_indices[i], values[i]);
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const NCPA::Vector<T> *values ) {
			this->assert_not_finalized();
			this->assert_row_in_range( rownum );
			this->assert_columns_match( values->size() );
			rows_[rownum]->copy( values );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) {
			this->assert_not_finalized();
			this->assert_column_in_range( colnum );
//			this->assert_row_in_range( nvals );
			this->clear_column( colnum );
			for (size_t i = 0; i < nvals; i++) {
				this->assert_row_in_range( row_indices[i] );
				rows_[row_indices[i]]->set( colnum, values[i] );
//				rows_[row_indices[i]][colnum] = values[i];
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<T> &values ) {
			this->assert_not_finalized();
			this->assert_column_in_range( colnum );
			this->assert_rows_match( values.size() );
			this->clear_column( colnum );
			for (size_t i = 0; i < rows(); i++) {
				if (values[i] != (T)0) {
					rows_[i]->set(colnum, values[i]);
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> &row_indices,
				const std::vector<T> &values ) {
			this->assert_not_finalized();
			this->assert_column_in_range( colnum );
			this->clear_column( colnum );
//			this->assert_row_in_range( row_indices.size() );
//			for (size_t i = 0; i < rows(); i++) {
//				auto cit = rows_[i].find( colnum );
//				if (cit != rows_[i].end()) {
//					rows_[i].erase(cit);
//				}
//			}
			for (size_t i = 0; i < values.size(); i++) {
				this->assert_row_in_range( row_indices[i] );
				rows_[row_indices[i]]->set(colnum, values[i]);
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const NCPA::Vector<T> *values ) {
			this->assert_not_finalized();
			this->assert_column_in_range( colnum );
			this->assert_rows_match( values->size() );
			this->clear_column( colnum );
			std::vector<size_t> inds = values->get_defined_indices();
			for (auto it = inds.cbegin(); it != inds.cend(); ++it) {
				this->assert_row_in_range( *it );
				rows_[*it]->set( colnum, values->get( *it ) );
			}
//
//
//			size_t this_ind = 0;
//			for (auto it = inds.cbegin(); it != inds.cend(); ++it) {
//				this->assert_row_in_range( *it );
//				// if the index is not defined in values, see if it is defined
//				// in this.  If it is, erase it.
//				while (this_ind < *it) {
//					auto cit = rows_[this_ind].find(colnum);
//					if (cit != rows_[this_ind].end()) {
//						rows_[this_ind].erase(cit);
//					}
//					this_ind++;
//				}
//				rows_[*it][colnum] = values->get( *it );
//			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* swap_rows( size_t r1, size_t r2 ) {
			this->assert_not_finalized();
			this->assert_row_in_range( r1 );
			this->assert_row_in_range( r2 );
			std::swap( this->rows_[r1], this->rows_[r2] );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* swap_columns( size_t c1, size_t c2 ) {
			this->assert_not_finalized();
			this->assert_column_in_range( c1 );
			this->assert_column_in_range( c2 );
			NCPA::Vector<T> *col1 = nullptr, *col2 = nullptr;
			this->get_column( c1, col1 );
			this->get_column( c2, col2 );
			this->set_column( c1, col2 )->set_column( c2, col1 );
			delete col1;
			delete col2;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

//		std::vector<T> solve( std::vector<T> b ) {
//			if (this->lu_lower_ == nullptr || this->lu_upper_ == nullptr) {
//				this->clear_lu_cache_();
//				this->LU_decompose();
//			}
//			if (!this->is_square()) {
//				throw std::domain_error(
//						"solve() not implemented for non-square matrices" );
//			}
//			std::vector<T> x;
//			NCPA::Matrix<T>::solve_with_LU_decomposition(
//					this->lu_lower_,
//					this->lu_upper_,
//					this->lu_permutation_,
//					b, x );
//
//			return x;
//		}

//		virtual void solve( NCPA::Vector<T> *b, NCPA::Vector<T> *solution ) {
//			std::vector<T> stdsol = this->solve( b->as_std() );
//			if (solution == nullptr) {
//				solution = NCPA::Vector<T>::build( NCPA::vector_t::SPARSE, stdsol );
//			} else {
//				solution->clear();
//				solution->from_std( stdsol );
//			}
//		}

		virtual bool is_ready() const {
			return (ncols_ > 0) && (rows_.size() > 0);
		}

		virtual NCPA::Matrix<T> *add( const NCPA::Matrix<T> *second ) {
			this->assert_not_finalized();
			this->assert_sizes_match( second );
			for (size_t i = 0; i < this->rows(); i++) {
				rows_[i]->add( second->get_row(i) );
//				NCPA::Vector<T> *row = nullptr;
//				second->get_row(i,row);
//				for (auto it = row->cbegin(); it != row->cend(); ++it) {
//					this->rows_[i][it->first] += it->second;
//				}
//				delete row;
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> *subtract( const NCPA::Matrix<T> *second ) {
			this->assert_not_finalized();
			this->assert_sizes_match( second );
			for (size_t i = 0; i < second->rows(); i++) {
				rows_[i]->subtract( second->get_row(i) );
//				NCPA::Vector<T> *row = nullptr;
//				second->get_row(i, row);
//				for (auto it = row->cbegin(); it != row->cend(); ++it) {
//					this->rows_[i][it->first] -= it->second;
//				}
//				delete row;
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> *transpose() {
			this->assert_not_finalized();
			NCPA::BasicMatrix<T> other( this->columns(), this->rows() );
			size_t counter = 0;
			for (auto rit = rows_.begin(); rit != rows_.end(); ++rit) {
				other.set_column(counter++, *rit);
			}
			swap( *this, other );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* multiply( const NCPA::Matrix<T> *second ) const {
			NCPA::Matrix<T> *product = new NCPA::BasicMatrix<T>();
			if (auto* derived = dynamic_cast<const NCPA::BasicMatrix<T>*>(second)) {
				NCPA::BasicMatrix<T>::multiply( this, derived, product );
			} else {
				NCPA::Matrix<T>::multiply( this, second, product );
			}
			return product;
		}


		static void multiply(
				const NCPA::BasicMatrix<T> *first,
				const NCPA::BasicMatrix<T> *second,
				NCPA::Matrix<T> *&product ) {
//			std::cout << "Calling multiply( BasicMatrix, BasicMatrix )" << std::endl;
			if (first->columns() != second->rows()) {
				std::ostringstream oss;
				oss << "Size mismatch in matrix product: Columns in first matrix ("
						<< first->columns() << ") must equal rows in second matrix ("
						<< second->rows() << ")";
				throw std::out_of_range( oss.str() );
			}
			if (product == nullptr) {
				product = new NCPA::BasicMatrix<T>( first->rows(), second->columns() );
			} else {
				product->initialize( first->rows(), second->columns() );
			}

			for (size_t r = 0; r < first->rows(); r++) {
				for (size_t c = 0; c < second->columns(); c++) {
//					NCPA::Vector<T> *row = nullptr, *col = nullptr;
//					first->get_row(r, row);
//					second->get_column(c, col);
					product->set( r, c,
							first->get_row(r)->scalar_product( second->get_column( c ) ) );
//					product->set( r, c, row->scalar_product( *col ) );
//					delete row;
//					delete col;
//					auto rit = row.cbegin();
//					auto cit = col.cbegin();
//					T prod = 0;
//					while (rit != row.cend() || cit != col.cend()) {
//						// product will only be nonzero if the indices match
//						if (rit->first == cit->first) {
//							prod += rit->second * cit->second;
//						} else if (rit->first > cit->first) {
//							++cit;
//						} else {
//							++rit;
//						}
//					}
//					product->set( r, c, prod );
				}
			}
		}

		virtual NCPA::Matrix<T> *scale( T factor ) {
			this->assert_not_finalized();
			for (auto rit = rows_.begin(); rit != rows_.end(); ++rit) {
				(*rit)->scale( factor );
//				for (auto cit = rit->begin(); cit != rit->end(); ++cit) {
//					cit->second *= factor;
//				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual std::vector<size_t> get_column_indices(size_t row) const {
			this->assert_row_in_range( row );
			return rows_[row]->get_defined_indices();
//			this->assert_finalized();
//			std::vector<size_t> inds;
//			for (auto it = rows_[row].cbegin(); it != rows_[row].cend(); ++it) {
//				inds.push_back( it->first );
//			}
//			std::sort( inds.begin(), inds.end() );
//			return inds;
		}

		virtual bool is_unequal( const NCPA::Matrix<T> &other ) const {
			if (this->rows() != other.rows()) {
				return true;
			}
			if (this->columns() != other.columns()) {
				return true;
			}
			for (size_t nrow = 0; nrow < rows(); nrow++) {
				if (rows_[nrow]->is_unequal( other.get_row( nrow ) ) ) {
					return true;
				}
//				std::vector<size_t> u_inds(2*this->columns()),
//						thiscols = this->get_column_indices(nrow),
//						othercols = other.get_column_indices(nrow);
//				std::vector<size_t>::iterator ind_it =
//						std::set_union(
//								thiscols.cbegin(),
//								thiscols.cend(),
//								othercols.cbegin(),
//								othercols.cend(),
//								u_inds.begin() );
//				u_inds.resize( ind_it - u_inds.begin() );
//				for (ind_it = u_inds.begin();
//						ind_it != u_inds.end();
//						++ind_it) {
//					if (this->get(nrow,*ind_it) != other.get(nrow,*ind_it)) {
//						return true;
//					}
//				}
			}
			return false;
		}

		// operators
		friend bool operator!=( const NCPA::BasicMatrix<T> &a, const NCPA::BasicMatrix<T> &b ) {
			return a.is_unequal( b );
		}
		friend bool operator==( const NCPA::BasicMatrix<T> &a, const NCPA::BasicMatrix<T> &b ) {
			return !(a.is_unequal( b ));
		}

	protected:

//		virtual void clear_lu_cache_() {
//			if (this->lu_lower_ != nullptr) {
//				delete this->lu_lower_;
//				this->lu_lower_ = nullptr;
//			}
//			if (this->lu_upper_ != nullptr) {
//				delete this->lu_upper_;
//				this->lu_upper_ = nullptr;
//			}
//			if (this->lu_permutation_ != nullptr) {
//				delete this->lu_permutation_;
//				this->lu_permutation_ = nullptr;
//			}
//		}

		std::vector<NCPA::BasicVector<T>*> rows_;
		mutable std::map<size_t,NCPA::BasicVector<T>*> cols_cache_;
		size_t ncols_ = 0;

		// cache variables
//		NCPA::DenseMatrix<T> *lu_lower_ = nullptr, *lu_upper_ = nullptr;
//		NCPA::DenseMatrix<int> *lu_permutation_ = nullptr;
	};
}

template<typename T>
void swap( NCPA::BasicMatrix<T> &a, NCPA::BasicMatrix<T> &b ) noexcept {
	using std::swap;
	::swap( static_cast<NCPA::Matrix<T> &>(a), static_cast<NCPA::Matrix<T> &>(b) );
	swap( a.rows_, b.rows_ );
	swap( a.ncols_, b.ncols_ );
	swap( a.cols_cache_, b.cols_cache_ );
//	swap( a.lu_lower_, b.lu_lower_ );
//	swap( a.lu_upper_, b.lu_upper_ );
//	swap( a.lu_permutation_, b.lu_permutation_ );
}


#endif
