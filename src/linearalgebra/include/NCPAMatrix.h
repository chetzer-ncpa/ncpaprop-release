/*
NCPAMatrix.h

This module provides the abstract base class template for Matrix objects for use in the NCPA
Linear Algebra library.

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
 * @file NCPAMatrix.h
 * @author Claus Hetzer
 * @version 1.0.0
 * @date 2023-11-10
 * @brief Provides the API for the NCPA::Matrix class.
 */

#ifndef NCPA_LINEARALGEBRA__MATRIX_H_INCLUDED
#define NCPA_LINEARALGEBRA__MATRIX_H_INCLUDED

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <type_traits>

#include "NCPACommon.h"
#include "NCPAVector.h"


namespace NCPA {
	/**
	 * An enumerated type for different kinds of matrices.
	 */
//	enum class matrix_t : unsigned int {
//		UNDEFINED = 0,		/**< Generic, undefined matrix */
//		DENSE,				/**< Standard dense (i.e. not storage-optimized) matrix */
//		SPARSE,				/**< Sparse, storage-optimized matrix */
//		DIAGONAL,			/**< Diagonal matrix */
//		TRIDIAGONAL			/**< Tridiagonal matrix */
//	};

	template<typename T> class Matrix;

}


/**
 * Swaps two Matrix objects.
 * @param a First object to swap
 * @param b Second object to swap
 */
template<typename T>
void swap( NCPA::Matrix<T> &a, NCPA::Matrix<T> &b ) noexcept;

template<typename T>
bool operator==(const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b);
template<typename T>
bool operator!=(const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b);

namespace NCPA {



	template<typename T>
	class Matrix {

	// prevent size_t and the like, matrix elements must be capable of being signed
	static_assert(!(std::is_unsigned<T>::value), "Matrix cannot be instantiated with unsigned arithmetic types.");

	public:
		Matrix() {}
		Matrix( const NCPA::Matrix<T> &other ) {
			this->finalized_ = other.finalized_;
		}
		Matrix( NCPA::Matrix<T> &&other ) noexcept {
			::swap( *this, other );
		}
		virtual ~Matrix() {}
		friend void ::swap<T>( NCPA::Matrix<T> &a, NCPA::Matrix<T> &b ) noexcept;


		/**
		 * @defgroup ncpa-matrix-api NCPA::Matrix API
		 * @{
		 */

		/**
		 * Clear contents and set dimensions of matrix.
		 * @param nrows The number of rows in the new matrix
		 * @param ncols The number of columns in the new matrix
		 * @return A pointer to the matrix
		 */
		virtual NCPA::Matrix<T> * initialize( size_t nrows, size_t ncols ) = 0;

		/**
		 * Creates a dynamically-allocated copy of this matrix.
		 * @return A pointer to the new matrix
		 */
		virtual NCPA::Matrix<T>* clone() const = 0;

		/**
		 * Clears any cached contents.
		 * @return A pointer to the matrix
		 */
		virtual NCPA::Matrix<T> * clear_cache() = 0;

		/**
		 * The number of rows in the matrix.
		 * @return The number of rows in the matrix.
		 */
		virtual size_t rows() const = 0;

		/**
		 * The number of columns in the matrix.
		 * @return The number of columns in the matrix.
		 */
		virtual size_t columns() const = 0;

		/**
		 * Returns an element of the matrix
		 * @param r The row number to get
		 * @param c The column number to get
		 * @return The contents of matrix[r,c]
		 */
		virtual T get( size_t r, size_t c ) const = 0;

		/**
		 * Sets an element of the matrix
		 * @param r The row number to set
		 * @param c The column number to set
		 * @param val The value to put at element [r,c]
		 * @return A pointer to the matrix
		 */
		virtual NCPA::Matrix<T>* set( size_t r, size_t c, T val ) = 0;

		/**
		 * Marks the matrix as finalized and performs any under-the-hood
		 * calculations that are needed for efficiency.
		 * @return A pointer to the matrix
		 */
		virtual NCPA::Matrix<T>* finalize() = 0;

		/**
		 * Zero out a row of the matrix.
		 * @param[in] rownum The number of the row to clear, indexed to 0.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* clear_row( size_t rownum ) = 0;

		/**
		 * Returns a row of the matrix.
		 * @param[in] rownum The number of the row to return, indexed to 0.
		 * @param[in,out] n The number of elements in the row.  If zero, will be set properly.  If nonzero,
		 * 		must match the value of <code>rows()</code> and the number of elements in <code>row</code>.
		 * @param[in,out] row The array in which to return the row values.  If <code>nullptr</code>,
		 * 		will be dynamically allocated.
		 */
		virtual void get_row( size_t rownum, size_t &n, T* &row ) const = 0;

		/**
		 * Returns a row of the matrix.
		 * @param[in] rownum The number of the row to return, indexed to 0.
		 * @param[out] row The vector in which to return the row values.  Any existing
		 * 		contents will be destroyed.
		 */
		virtual void get_row( size_t rownum, std::vector<T> &row ) const = 0;

		/**
		 * Returns a row of the matrix.
		 * @param[in] rownum The number of the row to return, indexed to 0.
		 * @param[out] row The vector in which to return the row values.  Any existing
		 * 		contents will be destroyed.
		 */
		virtual void get_row( size_t rownum, NCPA::Vector<T>* row ) const = 0;

		/**
		 * Returns a row of the matrix.
		 * @param[in] rownum The number of the row to return, indexed to 0.
		 * @return A vector containing the row values.
		 */
		virtual const NCPA::Vector<T>* get_row( size_t rownum ) const = 0;

		/**
		 * Zero out a column of the matrix.
		 * @param[in] colnum The number of the column to clear, indexed to 0.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* clear_column( size_t colnum ) = 0;

		/**
		 * Returns a column of the matrix.
		 * @param[in] colnum The number of the column to return, indexed to 0.
		 * @param[in,out] n The number of elements in the column.  If zero, will be set properly.  If nonzero,
		 * 		must match the value of <code>columns()</code> and the number of elements in <code>column</code>.
		 * @param[in,out] column The array in which to return the column values.  If <code>nullptr</code>,
		 * 		will be dynamically allocated.
		 */
		virtual void get_column( size_t colnum, size_t &n, T* &column ) const = 0;

		/**
		 * Returns a column of the matrix.
		 * @param[in] colnum The number of the column to return, indexed to 0.
		 * @param[out] column The vector in which to return the column values.  Any existing
		 * 		contents will be destroyed.
		 */
		virtual void get_column( size_t colnum, std::vector<T> &column ) const = 0;

		/**
		 * Returns a column of the matrix.
		 * This version is primarily of use working with instances of SparseMatrix, but
		 * could be useful in other circumstances as well.
		 * @param[in] colnum The number of the column to return, indexed to 0.
		 * @param[out] column The vector in which to return the column values.  Any existing
		 * 		contents will be destroyed.
		 */
		virtual void get_column( size_t colnum, NCPA::Vector<T> *column ) const = 0;

		/**
		 * Returns a column of the matrix.
		 * @param[in] colnum The number of the column to return, indexed to 0.
		 * @return A vector containing the column values.
		 */
		virtual const NCPA::Vector<T>* get_column( size_t colnum ) const = 0;

		/**
		 * Returns the defined column indices for a row.
		 * For the given row, returns a vector of the column indices for which a value has been
		 * defined.  That value may be zero.  This is primarily of use for storage-optimized
		 * implementations where zero values are not stored in memory.
		 */
		virtual std::vector<size_t> get_column_indices(size_t rownum) const = 0;

		/**
		 * Sets the values of specified columns in the specified row.
		 * @param[in] rownum The number of the row to set, indexed to 0.
		 * @param[in] nvals The number of values to set, which may be less than the row length.
		 * @param[in] column_indices The indices of the <code>nvals</code> columns to set.
		 * @param[in] values The values to set the <code>nvals</code> columns to.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_row( size_t rownum, size_t nvals,
				const size_t* column_indices, const T* values ) = 0;

		/**
		 * Sets the values of the specified row.
		 * @param[in] rownum The number of the row to set, indexed to 0.
		 * @param[in] values The values to set the row to.  Must contain <code>columns()</code> elements.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_row( size_t rownum, const std::vector<T> &values ) = 0;

		/**
		 * Sets the values of specified columns in the specified row.
		 * @param[in] rownum The number of the row to set, indexed to 0.
		 * @param[in] column_indices The indices of the columns to set.
		 * @param[in] values The values to set the columns to.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_row( size_t rownum,
				const std::vector<size_t> &column_indices,
				const std::vector<T> &values ) = 0;

		/**
		 * Sets the values of specified columns in the specified row.
		 * @param[in] rownum The number of the row to set, indexed to 0.
		 * @param[in] values The values to set the columns to.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_row( size_t rownum,
				const NCPA::Vector<T> *values ) = 0;

		/**
		 * Sets the values of specified rows in the specified column.
		 * @param[in] colnum The number of the column to set, indexed to 0.
		 * @param[in] nvals The number of values to set, which may be less than the column length.
		 * @param[in] row_indices The indices of the <code>nvals</code> rows to set.
		 * @param[in] values The values to set the <code>nvals</code> rows to.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) = 0;

		/**
		 * Sets the values of the specified column.
		 * @param[in] colnum The number of the column to set, indexed to 0.
		 * @param[in] values The values to set the column to.  Must contain <code>rows()</code> elements.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_column( size_t rownum,
				const std::vector<T> &values ) = 0;

		/**
		 * Sets the values of specified rows in the specified column.
		 * @param[in] colnum The number of the column to set, indexed to 0.
		 * @param[in] row_indices The indices of the rows to set.
		 * @param[in] values The values to set the rows to.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> &row_indices,
				const std::vector<T> &values ) = 0;

		/**
		 * Sets the values of specified rows in the specified column.
		 * @param[in] colnum The number of the column to set, indexed to 0.
		 * @param[in] values The values to set the rows to.
		 * @return A pointer to the matrix.
		 */
		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const NCPA::Vector<T> *values ) = 0;

		/**
		 * Swap two rows in the matrix.
		 * @param[in] r1 The first row number, indexed to 0.
		 * @param[in] r2 The second row number, indexed to 0.
		 */
		virtual NCPA::Matrix<T>* swap_rows( size_t r1, size_t r2 ) = 0;

		/**
		 * Swap two columns in the matrix.
		 * @param[in] c1 The first column number, indexed to 0.
		 * @param[in] c2 The second column number, indexed to 0.
		 */
		virtual NCPA::Matrix<T>* swap_columns( size_t c1, size_t c2 ) = 0;

		/**
		 * Tests inequality of two matrices in an implementation-specific method.
		 * @param other The matrix against which to test for inequality.
		 * @return <code>true</code> if any differences are found, <code>false</code> otherwise.
		 */
		virtual bool is_unequal( const NCPA::Matrix<T> &other ) const = 0;

		/**
		 * Adds another matrix to this matrix.
		 * @param[in] second The matrix to add to this one.  Matrix dimensions must agree.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T> *add( const NCPA::Matrix<T> *second ) = 0;

		/**
		 * Subtracts another matrix from this matrix.
		 * @param[in] second The matrix to subtract from this one.  Matrix dimensions must agree.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T> *subtract( const NCPA::Matrix<T> *second ) = 0;

		/**
		 * Multiplies another matrix by this matrix.
		 * @param[in] second The RHS matrix.
		 * @param[out] product The product, which will be dynamically allocated
		 * 		if it is <code>nullptr</code>.
		 */
		virtual NCPA::Matrix<T> * multiply( const NCPA::Matrix<T> *second ) const = 0;

		/**
		 * Transposes this matrix.
		 * @return A pointer to this matrix after transposition.
		 */
		virtual NCPA::Matrix<T> *transpose() = 0;

		/**
		 * Multiplies this matrix by a scalar value.
		 * @param[in] factor The scalar by which to multiply all elements of the matrix.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T> *scale( T factor ) = 0;

		/**
		 * Perform an LU decomposition, caching the results.
		 * @return A pointer to this matrix.
		 */
//		virtual NCPA::Matrix<T> *LU_decompose() = 0;

		/**
		 * Return the results of an LU decomposition.
		 * @param[out] lower The lower-triangular portion of the decomposition.
		 * @param[out] upper The upper-triangular portion of the decomposition.
		 * @param[out] permut The permutation matrix such that permut*this == lower*upper
		 */
//		virtual void get_LU_decomposition(
//				NCPA::Matrix<T> *&lower,
//				NCPA::Matrix<T> *&upper,
//				NCPA::Matrix<int> *&permut ) const = 0;

//		/**
//		 * Solve the matrix equation <code><b>this</b>*x = b</code> for x.
//		 * The details of the solution are implementation-specific.
//		 * @param[in] b The right-hand side of the matrix equation.
//		 * @return The x vector in <code><b>this</b>*x = b</code>.
//		 */
//		virtual std::vector<T> solve( std::vector<T> b ) = 0;
//
//		/**
//		 * Solve the matrix equation <code><b>this</b>*x = b</code> for x.
//		 * The details of the solution are implementation-specific.
//		 * @param[in] b The right-hand side of the matrix equation.
//		 * @return The x vector in <code><b>this</b>*x = b</code>.
//		 */
//		virtual void solve( NCPA::Vector<T> *b, NCPA::Vector<T> *solution ) = 0;



		/**
		 * @}
		 */

		// ------------------------------------------------------------------------------------
		// concrete base methods below here.  Base classes may implement their own
		// versions
		// ------------------------------------------------------------------------------------

		/**
		 * Asserts that the matrix has been finalized and can be used for
		 * computations, as determined by the implementation.
		 * @throws std::logic_error if the matrix is not finalized
		 */
		virtual void assert_finalized() const {
			if (!this->is_finalized()) {
				throw std::logic_error( "Matrix not finalized!" );
			}
		}

		/**
		 * Asserts that the matrix has not been finalized and can be modified.
		 * @throws std::logic_error if the matrix is finalized
		 */
		virtual void assert_not_finalized() const {
			if (this->is_finalized()) {
				throw std::logic_error( "Matrix already finalized!" );
			}
		}

		/**
		 * Asserts that the matrix is square.
		 * @throws std::invalid_argument if the matrix is not square
		 */
		virtual void assert_square() const {
			if (!this->is_square()) {
				throw std::invalid_argument( "Matrix is not square!" );
			}
		}

		/**
		 * Asserts that the matrix has a certain number of rows.
		 * @param r The number of rows the matrix should have.
		 * @throws std::out_of_range if the matrix does not have that many rows.
		 */
		virtual void assert_rows_match( size_t r ) const {
			if (r != this->rows()) {
				std::ostringstream oss;
				oss << "Size mismatch: Matrix has " << this->rows()
						<< " rows, but operation has " << r << " elements.";
				throw std::out_of_range( oss.str() );
			}
		}

		/**
		 * Asserts that the matrix has a certain number of columns.
		 * @param c The number of columns the matrix should have.
		 * @throws std::out_of_range if the matrix does not have that many columns.
		 */
		virtual void assert_columns_match( size_t c ) const {
			if (c != this->columns()) {
				std::ostringstream oss;
				oss << "Size mismatch: Matrix has " << this->columns()
						<< " columns, but operation has " << c << " elements.";
				throw std::out_of_range( oss.str() );
			}
		}

		/**
		 * Asserts that the matrix has a certain number of rows and columns.
		 * @param r The number of rows the matrix should have.
		 * @param c The number of columns the matrix should have.
		 * @throws std::out_of_range if the matrix dimensions do not match.
		 */
		virtual void assert_sizes_match( size_t r, size_t c ) const {
			this->assert_rows_match( r );
			this->assert_columns_match( c );
		}

		/**
		 * Asserts that a row number is valid for this matrix.
		 * @param r The index of the row (starting at zero).
		 * @throws std::out_of_range if the matrix does not have that many rows.
		 */
		virtual void assert_row_in_range( size_t r ) const {
			if (r >= this->rows()) {
				std::ostringstream oss;
				oss << "Requested row " << r
						<< " is out of range for Matrix with "
						<< rows() << " rows.";
				throw std::out_of_range( oss.str() );
			}
		}

		/**
		 * Asserts that a column number is valid for this matrix.
		 * @param c The index of the column (starting at zero).
		 * @throws std::out_of_range if the matrix does not have that many columns.
		 */
		virtual void assert_column_in_range( size_t c ) const {
			if (c >= this->columns()) {
				std::ostringstream oss;
				oss << "Requested column " << c
						<< " is out of range for Matrix with "
						<< columns() << " columns.";
				throw std::out_of_range( oss.str() );
			}
		}

		/**
		 * Asserts that an element index is valid for this matrix.
		 * @param r The index of the row (starting at zero).
		 * @param c The index of the column (starting at zero).
		 * @throws std::out_of_range if the matrix does not have that many rows or columns.
		 */
		virtual void assert_element_in_range( size_t r, size_t c ) const {
			this->assert_row_in_range( r );
			this->assert_column_in_range( c );
		}

		/**
		 * Asserts that another matrix has the same dimensions.
		 * @param other The other matrix to test against.
		 * @throws std::out_of_range if the dimensions of the two matrices do not match.
		 */
		virtual void assert_sizes_match( const NCPA::Matrix<T> *other ) const {
			if (this->rows() != other->rows() || this->columns() != other->columns()) {
				std::ostringstream oss;
				oss << "Matrix size mismatch between RHS matrix of size "
						<< other->rows() << "x" << other->columns()
						<< " and LHS matrix of size " << this->rows() << "x"
						<< this->columns() << std::endl;
				throw std::out_of_range( oss.str() );
			}
		}

//		/**
//		 * Lock the matrix against further modifications.
//		 * Sets the matrix to locked, so that calling the <code>assert_not_finalized()</code> method
//		 * will throw an exception.
//		 * @return A pointer to this matrix.
//		 */
//		virtual NCPA::Matrix<T> * lock() {
//			locked_ = true;
//			return this;
//		}
//
//		/**
//		 * Unlock the matrix for further modifications.
//		 * Sets the matrix to unlocked, so that calling the <code>assert_not_finalized()</code> method
//		 * will throw an exception, and clears any cached data.
//		 * @return A pointer to this matrix.
//		 */
//		virtual NCPA::Matrix<T> * unlock() {
//			this->clear_cache();
//			locked_ = false;
//			return this;
//		}
//
//		/**
//		 * Returns whether the matrix is locked or not.
//		 * @return <code>true</code> if the matrix is locked, <code>false</code> otherwise.
//		 */
//		virtual bool locked() const {
//			return locked_;
//		}
//
//
//		/**
//		 * Assert that the matrix is locked, unless locking is permissive.
//		 * If unlocked and locking is strict (i.e. <code>is_lock_enforcing()</code> returns <code>true</code>,
//		 * throws std::logic_error.  If unlocked and locking is permissive, locks the matrix and returns.
//		 * If locked, returns.
//		 */
//		virtual void assert_locked() {
//			if (!this->locked()) {
//				if (this->is_lock_enforcing()) {
//					throw std::logic_error( "Matrix must be locked to perform this operation!" );
//				} else {
//					this->lock();
//				}
//			}
//		}
//
//		/**
//		 * Assert that the matrix is unlocked, unless locking is permissive.
//		 * If locked and locking is strict (i.e. <code>is_lock_enforcing()</code> returns <code>true</code>,
//		 * throws std::logic_error.  If locked and locking is permissive, unlocks the matrix and returns.
//		 * If locked, returns.
//		 */
//		virtual void assert_not_finalized() {
//			if (this->is_lock_enforcing()) {
//				if (this->locked()) {
//					throw std::logic_error( "Matrix is locked!" );
//				}
//			} else {
//				if (this->locked()) {
//					this->unlock();
//				}
//			}
//		}
//
//
//		/**
//		 * Toggles strict enforcement of the lock.
//		 * By default locking is strict, meaning that an operation that requires a locked matrix, when it
//		 * calls the <code>enforce_lock()</code> method, will throw an exception if the matrix has not
//		 * been explicitly locked (i.e. the user has acknowledged that cached values may depend on the
//		 * current state of the matrix and it should not be changed after this point).  This function
//		 * can set enforcement to permissive, meaning that calling <code>enforce_lock()</code> on an
//		 * unlocked matrix will lock it and proceed from there.
//		 * @param tf The value to set for strict locking (i.e. <code>true</code> means strict,
//		 * 		<code>false</code> means permissive.
//		 * @return A pointer to this matrix.
//		 */
//		virtual NCPA::Matrix<T> * enforce_lock( bool tf ) {
//			this->enforce_lock_ = tf;
//			return this;
//		}
//
//		/**
//		 * Returns whether the matrix is set to strict lock enforcement.
//		 * @return <code>true</code> if strict enforcement is set, <code>false</code> otherwise.
//		 * @see NCPA::Matrix<T>::enforce_lock()
//		 */
//		virtual bool is_lock_enforcing() const {
//			return this->enforce_lock_;
//		}

		/**
		 * Calculates an LU decomposition using Gaussian elimination with optional partial pivoting.
		 * Subclasses may override this method if appropriate to the implementation.
		 * @param[in] base The matrix to decompose.
		 * @param[out] lower The lower-triangular part of the decomposition with diagonal=1.
		 * @param[out] upper The upper-triangular part of the decomposition.
		 * @param[out] permutation The permutation matrix such that <code>permutation*base == lower*upper</code>.
		 * @param[in] tol The tolerance to use to determine degeneracy.
		 */
//		static void calculate_LU_decomposition(
//				NCPA::Matrix<T> *base, NCPA::Matrix<T> *lower, NCPA::Matrix<T> *upper,
//				NCPA::Matrix<int> *permutation, bool pivot = true,
//				T tol = ((T)1e-20) ) {
//			base->assert_ready();
//			base->assert_square();
//			base->assert_locked();
//
//			size_t i, j, k;
//			size_t N = base->rows();
//			T maxA, absA;
//
//			upper->initialize( N, N );
//			lower->initialize( N, N );
//			permutation->identity( N, N );
//			upper->copy( base );
//
//			for (k = 0; k < N; k++) {
//				if (pivot) {
//					size_t pivotRow = k;
//					T maxVal = upper->get( k, k );
//
//					// find pivot row and swap
//					for (i = k+1; i < N; i++) {
//						T x = upper->get( i, k );
//						if (std::abs(x) > std::abs(maxVal)) {
//							maxVal = x;
//							pivotRow = i;
//						}
//					}
//					if (std::abs(maxVal) < std::abs(tol)) {
//						throw std::runtime_error( "Matrix is degenerate" );
//					}
//					upper->swap_rows( k, pivotRow );
//					permutation->swap_rows( k, pivotRow );
//					lower->swap_rows( k, pivotRow );
//				}
//
//				// perform elimination
//				for (i = k+1; i < N; i++) {
//					T ratio = upper->get( i, k ) / upper->get( k, k );
//					lower->set( i, k, ratio );
//					for (j = 0; j < N; j++) {
//						T diff = upper->get( i, j ) - ratio * upper->get( k, j );
//						upper->set( i, j, diff );
//					}
//				}
//			}
//
//			for (k = 0; k < N; k++) {
//				lower->set( k, k, ((T)1.0) );
//			}
//		}

		/**
		 * Casts the elements of this matrix to elements of a matrix of a different type.
		 * @param[in,out] to The matrix to which the elements of this matrix should be cast.  The target
		 * 		Matrix will be initialized before casting.
		 */
		template<typename U>
		void convert( NCPA::Matrix<U> *to ) {
			to->initialize( this->rows(), this->columns() );
			for (size_t r = 0; r < this->rows(); r++) {
				for (size_t c = 0; c < this->columns(); c++) {
					to->set( r, c, ((U)(this->get(r,c))) );
				}
			}
		}

		/**
		 * Copies the elements of another Matrix of the same type to this Matrix.
		 * Essentially this is meant for the reuse of existing matrices without explicitly going through
		 * the dynamic memory deallocation/allocation process.
		 * @param[in] source The matrix to copy element-wise.
		 * @return A pointer to this matrix with the new elements.
		 */
		virtual NCPA::Matrix<T> *copy( const NCPA::Matrix<T> *source ) {
			this->initialize( source->rows(), source->columns() );
			for (size_t r = 0; r < source->rows(); r++) {
				for (size_t c = 0; c < source->columns(); c++) {
					this->set( r, c, source->get( r, c ) );
				}
			}
			return this;
		}

		/**
		 * Returns the diagonal of this matrix.
		 * @return The values along the diagonal of the matrix.
		 */
		virtual std::vector<T> get_diagonal() const {
			size_t n = rows() < columns() ? rows() : columns();
			std::vector<T> diag(n);
			for (size_t ii = 0; ii < n; ii++) {
				diag[ii] = this->get(ii,ii);
			}
			return diag;
		}

		/**
		 * Returns the diagonal of this matrix.
		 * @param[out] The number of values in the diagonal.  Will be set if 0, otherwise must match the
		 * 		actual length of the diagonal.
		 * @param[out] The values of the diagonal.  Will be dynamically allocated if <code>nullptr</code>.
		 * @throws std::out_of_range if the size is nonzero and incorrect.
		 */
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

		/**
		 * Returns an off-diagonal of this matrix.
		 * Returns the m'th off-diagonal, where m=abs(offset), taken from the lower half if
		 * offset is negative, and the upper half if offset is positive.
		 * @param[in] offset The number of the off-diagonal to return.  Negative <code>offset</code>
		 * 		returns from the lower half, positive from the upper half.
		 * @param[in,out] n The number of elements returned, which is size(diagonal)-abs(offset).  Will
		 * 		be set if zero, otherwise will throw an exception if incorrect.
		 * @param[out] values The array in which to store the values.
		 * @throws std::out_of_range if <code>n</code> is nonzero and incorrect.
		 */
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

		/**
		 * Returns an off-diagonal of this matrix.
		 * Returns the m'th off-diagonal, where m=abs(offset), taken from the lower half if
		 * offset is negative, and the upper half if offset is positive.
		 * @param[in] offset The number of the off-diagonal to return.  Negative <code>offset</code>
		 * 		returns from the lower half, positive from the upper half.
		 * @return A vector containing the off-diagonal requested.
		 */
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

		/**
		 * Reconfigures the matrix as an identity matrix (i.e. 1 on the diagonal, 0 elsewhere).
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T> * identity( size_t r, size_t c ) {
			this->assert_not_finalized();
			this->initialize( r, c );
			size_t N = NCPA::min<size_t>( r, c );
			for (size_t i = 0; i < N; i++) {
				this->set( i, i, ((T)1.0) );
			}
			return this;
		}

		/**
		 * Returns <code>true</code> if the matrix is diagonal.
		 * Uses a tolerance that defaults to 0 to identify zero elements.
		 * @param tol The tolerance to use (element is 0 if <= <code>tol</code>.
		 * @return <code>true</code> if the matrix is diagonal, <code>false</code> otherwise.
		 */
		bool is_diagonal(T tol = ((T)0)) const {
			for (size_t r = 0; r < this->rows(); r++) {
				for (size_t c = 0; c < this->columns(); c++) {
					if (r != c && std::abs(this->get(r,c)) > std::abs(tol)) {
						return false;
					}
				}
			}
			return true;
		}

		/**
		 * Returns whether the matrix has been finalized.
		 * @returns <code>true</code> if the matrix has been finalized, <code>false</code>
		 * 		otherwise.
		 */
		virtual bool is_finalized() const {
			return finalized_;
		}

		/**
		 * Returns <code>true</code> if the matrix is lower-triangular.
		 * Uses a tolerance that defaults to 0 to identify zero elements.
		 * @param tol The tolerance to use (element is 0 if <= <code>tol</code>.
		 * @return <code>true</code> if the matrix is lower-triangular, <code>false</code> otherwise.
		 */
		bool is_lower_triangular(T tol = ((T)0)) const {
			for (size_t r = 0; r < this->rows(); r++) {
				for (size_t c = r+1; c < this->columns(); c++) {
					if (std::abs(this->get(r,c)) > std::abs(tol)) {
						return false;
					}
				}
			}
			return true;
		}

		/**
		 * Returns whether the matrix is ready to be used.
		 * @returns <code>true</code> if the matrix can be used for operations, <code>false</code>
		 * 		otherwise.
		 */
//		virtual bool is_ready() const {
//			return ready_;
//		}

		/**
		 * Returns <code>true</code> if the matrix is square (<code>rows() == columns()</code>).
		 * @return <code>true</code> if the matrix is square, <code>false</code> otherwise.
		 */
		bool is_square() const {
			return (rows() == columns());
		}

		/**
		 * Returns <code>true</code> if the matrix is tridiagonal.
		 * Uses a tolerance that defaults to 0 to identify zero elements.
		 * @param tol The tolerance to use (element is 0 if <= <code>tol</code>.
		 * @return <code>true</code> if the matrix is tridiagonal, <code>false</code> otherwise.
		 */
		bool is_tridiagonal(T tol = ((T)0)) const {
			for (size_t r = 0; r < this->rows(); r++) {
				for (size_t c = 0; c < this->columns(); c++) {
					if (r < (c-1) && r > (c+1) && std::abs(this->get(r,c)) > std::abs(tol)) {
						return false;
					}
				}
			}
			return true;
		}

		/**
		 * Returns <code>true</code> if the matrix is upper-triangular.
		 * Uses a tolerance that defaults to 0 to identify zero elements.
		 * @param tol The tolerance to use (element is 0 if <= <code>tol</code>.
		 * @return <code>true</code> if the matrix is upper-triangular, <code>false</code> otherwise.
		 */
		bool is_upper_triangular(T tol = ((T)0)) const {
			for (size_t r = 0; r < this->rows(); r++) {
				for (size_t c = 0; c < r; c++) {
					if (std::abs(this->get(r,c)) > std::abs(tol)) {
						return false;
					}
				}
			}
			return true;
		}

		/**
		 * Multiply this matrix by another matrix and keep the product as this matrix.
		 * @param The matrix to multiply by.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* multiply( const NCPA::Matrix<T> *second,
				NCPA::Matrix<T> *&product ) const {
			NCPA::Matrix<T>::multiply( this, second, product );
			return product;
		}

		/**
		 * Multiply two matrices together.
		 * @param[in] first The left matrix.
		 * @param[in] second The right matrix.
		 * @param[out] product The new matrix generated.
		 * @throws std::invalid_argument if the matrix sizes are not compatible.
		 * @throws std::logic_error if a null pointer is provided for the product matrix.
		 */
		static void multiply(
				const NCPA::Matrix<T> *first,
				const NCPA::Matrix<T> *second,
				NCPA::Matrix<T> *&product ) {
			first->assert_finalized();
			second->assert_finalized();
			if (first->columns() != second->rows()) {
				std::ostringstream oss;
				oss << "Size mismatch in matrix product: Columns in first matrix ("
						<< first->columns() << ") must equal rows in second matrix ("
						<< second->rows() << ")";
				throw std::invalid_argument( oss.str() );
			}
			if (product == nullptr) {
				throw std::logic_error( "No product matrix provided" );
			}
			product->initialize(first->rows(), second->columns());
			std::vector<NCPA::Vector<T>*> secondcols( second->columns(), nullptr );
			for (size_t r = 0; r < first->rows(); r++) {
//				if (firstrows[r] == nullptr) {
//					firstrows[r] = first->get_row( r );
//				}
//				NCPA::Vector<T> *firstrow = first->get_row( r );
				for (size_t c = 0; c < second->columns(); c++) {
					product->set( r, c,
							first->get_row( r )->scalar_product(
									second->get_column( c ) ) );
//					if (secondcols[c] == nullptr) {
//						secondcols[c] = second->get_column( c );
//					}
//					product->set( r, c, firstrow->scalar_product(*secondcols[c]) );
//					if (r == first->rows()-1) {
//						delete secondcols[c];
//					}
				}
//				delete firstrow;
			}
		}

		/**
		 * Negates the matrix.
		 * Scales the matrix by a factor of -1.
		 * @return A pointer to this matrix.
		 * @throws std::logic_error if the type of the matrix is not signed.
		 */
		virtual NCPA::Matrix<T> *negative() {
			this->assert_not_finalized();
			return this->scale(-1.0);
		}

		/**
		 * Prints a text version of the matrix.
		 * @param o The output stream to print to.
		 */
		virtual void print( std::ostream &o = std::cout ) const {
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

		/**
		 * Sets the diagonal elements of the matrix.
		 * Can set the first <code>n</code> diagonal elements, where
		 * <code>n <= min( rows(), columns() )</code>.
		 * @param n The number of elements to set.
		 * @param values The values to emplace in the diagonal.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_diagonal( size_t n, const T* values ) {
			this->assert_not_finalized();
			this->assert_element_in_range( n-1, n-1 );
			for (size_t i = 0; i < n; i++) {
				this->set( i, i, values[i] );
			}
			return this;
		}

		/**
		 * Sets the diagonal elements of the matrix.
		 * Can set the first <code>n</code> diagonal elements, where
		 * <code>n <= min( rows(), columns() )</code>.
		 * @param values The values to emplace in the diagonal.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_diagonal( const std::vector<T> &values ) {
			this->assert_not_finalized();
			this->assert_element_in_range( values.size()-1, values.size()-1 );
			for (size_t i = 0; i < values.size(); i++) {
				this->set( i, i, values[i] );
			}
			return this;
		}

		/**
		 * Sets the diagonal elements of the matrix to a constant.
		 * @param val The value to emplace in the diagonal.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_diagonal( T val ) {
			this->assert_not_finalized();
			size_t n = this->calculate_diagonal_length_();
			for (size_t i = 0; i < n; i++) {
				this->set( i, i, val );
			}
			return this;
		}

		/**
		 * Sets one of the off-diagonals of the matrix to a constant.
		 * Sets the m'th off-diagonal, where m=abs(offset), taken from the lower half if
		 * offset is negative, and the upper half if offset is positive.
		 * @param[in] offset The number of the off-diagonal to set.  Negative <code>offset</code>
		 * 		sets in the lower half, positive in the upper half.
		 * @param val The value to emplace in the off-diagonal.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_offdiagonal( int offset, T val ) {
			this->assert_not_finalized();
			int ndiag;
			size_t rstart, cstart;
			calculate_offdiagonal_parameters_( offset, ndiag, rstart, cstart );
			for (size_t i = 0; i < ndiag; i++) {
				this->set( rstart+i, cstart+i, val );
			}
			return this;
		}

		/**
		 * Sets the values of one of the off-diagonals of the matrix.
		 * Sets the m'th off-diagonal, where m=abs(offset), taken from the lower half if
		 * offset is negative, and the upper half if offset is positive.
		 * Can set the first <code>n</code> diagonal elements, where
		 * <code>n <= min( rows(), columns() )</code>.
		 * @param[in] offset The number of the off-diagonal to set.  Negative <code>offset</code>
		 * 		sets in the lower half, positive in the upper half.
		 * @param n The number of values to set
		 * @param values The values to emplace in the off-diagonal.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_offdiagonal( int offset,
				size_t n, const T* values ) {
			this->assert_not_finalized();
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

		/**
		 * Sets the values of one of the off-diagonals of the matrix.
		 * Sets the m'th off-diagonal, where m=abs(offset), taken from the lower half if
		 * offset is negative, and the upper half if offset is positive.
		 * Can set the first <code>n</code> diagonal elements, where
		 * <code>n <= min( rows(), columns() )</code>.
		 * @param[in] offset The number of the off-diagonal to set.  Negative <code>offset</code>
		 * 		sets in the lower half, positive in the upper half.
		 * @param values The values to emplace in the off-diagonal.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_offdiagonal( int offset,
				const std::vector<T> &values ) {
			this->assert_not_finalized();
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

		/**
		 * Sets the values of specified columns in the specified row.
		 * @param[in] rownum The number of the row to set, indexed to 0.
		 * @param[in] column_indices The indices of the columns to set.
		 * @param[in] values The values to set the columns to.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						std::initializer_list<size_t> column_indices,
						std::initializer_list<T> values ) {
			return this->set_row( rownum,
					std::vector<size_t>(column_indices),
					std::vector<T>( values ) );
		}

		/**
		 * Sets the values of specified rows in the specified column.
		 * @param[in] colnum The number of the column to set, indexed to 0.
		 * @param[in] row_indices The indices of the rows to set.
		 * @param[in] values The values to set the rows to.
		 * @return A pointer to this matrix.
		 */
		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				std::initializer_list<size_t> row_indices,
				std::initializer_list<T> values ) {
			return this->set_column( colnum,
					std::vector<size_t>(row_indices),
					std::vector<T>( values ) );
		}


		/**
		 * Tests inequality of two Matrix objects.
		 * Uses the <code>is_unequal()</code> method to test for inequality.
		 * @param a The left matrix.
		 * @param b The right matrix.
		 * @return <code>true</code> if the matrix elements test equal, <code>false</code> otherwise.
		 */
		friend bool operator!=( const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b ) {
			return a.is_unequal( b );
		}

		/**
		 * Tests equality of two Matrix objects.
		 * Uses the <code>is_unequal()</code> method to test for equality.
		 * @param a The left matrix.
		 * @param b The right matrix.
		 * @return <code>true</code> if the matrix elements test unequal, <code>false</code> otherwise.
		 */
		friend bool operator==( const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &b ) {
			return !(a.is_unequal( b ));
		}

	protected:
//		bool locked_ = false;
//		bool enforce_lock_ = true;
		bool finalized_ = false;

		/**
		 * Calculates the length of the matrix diagonal.
		 * Returns <code>min( rows(), columns() )</code> to account for non-square matrices.
		 * @return The number of elements in the matrix diagonal.
		 */
		virtual size_t calculate_diagonal_length_() const {
			return NCPA::min<size_t>( this->rows(), this->columns() );
		}

		/**
		 * Calculates the length and start element of the specified off-diagonal.
		 * Calculates for the m'th off-diagonal, where m=abs(offset), taken from the lower half if
		 * offset is negative, and the upper half if offset is positive.
		 * @param[in] offset The number of the off-diagonal to set.  Negative <code>offset</code>
		 * 		sets in the lower half, positive in the upper half.
		 * @param[out] ndiag The number of points in the off-diagonal.
		 * @param[out] rstart The row number of the first off-diagonal element.
		 * @param[out] cstart The column number of the first off-diagonal element.
		 * @throws std::out_of_range if the offset is too large (i.e. the off-diagonal does not exist).
		 */
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


template<typename T>
void swap( NCPA::Matrix<T> &a, NCPA::Matrix<T> &b ) noexcept {
	using std::swap;
	swap( a.finalized_, b.finalized_ );
}



#endif
