#ifndef NCPA__LINEARALGEBRA_NCPATRIDIAGONALMATRIX_H_INCLUDED_
#define NCPA__LINEARALGEBRA_NCPATRIDIAGONALMATRIX_H_INCLUDED_

#include "NCPAMatrix.h"
#include "NCPASparseMatrix.h"
#include "NCPASparseVector.h"
#include "NCPACommon.h"
#include <vector>

namespace NCPA { template<typename T> class TridiagonalMatrix; }

template<typename T>
void swap( NCPA::TridiagonalMatrix<T> &a, NCPA::TridiagonalMatrix<T> &b ) noexcept;
template<typename T>
bool operator==(const NCPA::TridiagonalMatrix<T> &a, const NCPA::TridiagonalMatrix<T> &b);
template<typename T>
bool operator!=(const NCPA::TridiagonalMatrix<T> &a, const NCPA::TridiagonalMatrix<T> &b);

namespace NCPA {

	template<typename T>
	class TridiagonalMatrix : public NCPA::SparseMatrix<T> {
	public:
		TridiagonalMatrix() : NCPA::SparseMatrix<T>() {}

		TridiagonalMatrix( size_t nrows, size_t ncolumns )
				: NCPA::SparseMatrix<T>(nrows,ncolumns) {}

		TridiagonalMatrix( const TridiagonalMatrix<T> &other )
				: NCPA::SparseMatrix<T>( other ) {}

		TridiagonalMatrix( TridiagonalMatrix<T> &&other ) noexcept
				: NCPA::SparseMatrix<T>() {
			::swap( *this, other );
		}

		friend void ::swap<T>(
				NCPA::TridiagonalMatrix<T> &a, NCPA::TridiagonalMatrix<T> &b ) noexcept;

		TridiagonalMatrix<T>& operator=( NCPA::TridiagonalMatrix<T> other ) {
			::swap( *this, other );
			return *this;
		}

		virtual NCPA::Matrix<T> * clone() const {
			return static_cast<NCPA::Matrix<T> *>( new NCPA::SparseMatrix<T>( *this ) );
		}

		virtual NCPA::Matrix<T>* set( size_t row, size_t col, T val ) {
			this->check_off_diagonal_( row, col );
			return this->SparseMatrix<T>::set( row, col, val );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum, size_t nvals,
						const size_t * column_indices, const T* values ) {
			for (size_t i = 0; i < nvals; i++) {
				this->check_off_diagonal_( rownum, column_indices[i] );
			}
			return this->SparseMatrix<T>::set_row(
					rownum, nvals, column_indices, values );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<T> &values ) {
			size_t ilower, iupper;
			this->get_column_range_( rownum, ilower, iupper );
			if (values.size() <= 3) {
				NCPA::SparseVector<T> sv;
				for (size_t i = 0; i < values.size(); i++) {
					sv[i+ilower] = values[i];
				}
				return this->SparseMatrix<T>::set_row( rownum, sv );
			} else {
				std::vector<T> v( values.size() );
				for (size_t i = ilower; i <= iupper; i++) {
					v[i] = values[i];
				}
				return this->SparseMatrix<T>::set_row( rownum, v );
			}
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<size_t> &column_indices,
						const std::vector<T> &values ) {
			for (auto it = column_indices.cbegin(); it != column_indices.cend(); ++it) {
				this->check_off_diagonal_( rownum, *it );
			}
			return this->SparseMatrix<T>::set_row( rownum, column_indices, values );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const NCPA::SparseVector<T> &values ) {
			for (auto it = values.cbegin(); it != values.cend(); ++it) {
				this->check_off_diagonal_( rownum, it->first );
			}
			return this->SparseMatrix<T>::set_row(
					rownum, values );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) {
			for (size_t i = 0; i < nvals; i++) {
				this->check_off_diagonal_( row_indices[i], colnum );
			}
			return this->SparseMatrix<T>::set_column(
				colnum, nvals, row_indices, values );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<T> &values ) {
			size_t ilower, iupper;
			this->get_row_range_( colnum, ilower, iupper );
			if (values.size() <= 3) {
				NCPA::SparseVector<T> sv;
				for (size_t i = ilower; i <= iupper; i++) {
					sv[i] = values[i-ilower];
				}
				return this->SparseMatrix<T>::set_column(
						colnum, sv );
			} else {
				std::vector<T> v( values.size() );
				for (size_t i = ilower; i <= iupper; i++) {
					v[i] = values[i];
				}
				return this->SparseMatrix<T>::set_column(
						colnum, v );
			}
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> &row_indices,
				const std::vector<T> &values ) {
			for (auto it = row_indices.cbegin(); it != row_indices.cend(); ++it) {
				this->check_off_diagonal_( *it, colnum );
			}
			return this->SparseMatrix<T>::set_column(
					colnum, row_indices, values );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const NCPA::SparseVector<T> &values ) {
			for (auto it = values.cbegin(); it != values.cend(); ++it) {
				this->check_off_diagonal_( it->first, colnum );
			}
			return this->SparseMatrix<T>::set_column(
					colnum, values );
		}

		static void multiply(
				const NCPA::TridiagonalMatrix<T> *first,
				const NCPA::TridiagonalMatrix<T> *second,
				NCPA::Matrix<T> *&product ) {

			if (first->columns() != second->rows()) {
				std::ostringstream oss;
				oss << "Size mismatch in matrix product: Columns in first matrix ("
						<< first->columns() << ") must equal rows in second matrix ("
						<< second->rows() << ")";
				throw std::out_of_range( oss.str() );
			}
			product->initialize( first->rows(), second->columns() );

			size_t rmin, rmax, cmin, cmax, mindiag, maxdiag;
			for (size_t r = 0; r < first->rows(); r++) {
				mindiag = (size_t)NCPA::max<int>( (int)r - 2, 0 );
				maxdiag = (size_t)NCPA::min<int>( (int)r + 2, first->rows()-1 );
				for (size_t c = mindiag; c <= maxdiag; c++) {
					first->get_row_range_( c, rmin, rmax );
					T val = 0.0;
					for (size_t i = rmin; i <= rmax; i++) {
						val += first->get(r,i) * second->get(i,c);
					}
					product->set( r, c, val );
				}
			}
		}

		static NCPA::Matrix<T>* multiply(
				const NCPA::TridiagonalMatrix<T> *first,
				const NCPA::TridiagonalMatrix<T> *second ) {
			NCPA::Matrix<T> *product = new NCPA::SparseMatrix<T>();
			NCPA::TridiagonalMatrix<int>::multiply( first, second, product );
			return product;
		}

		virtual std::vector<size_t> get_column_indices(size_t row) const {
			this->assert_ready_();
			size_t cmin, cmax;
			this->get_column_range_( row, cmin, cmax );
			std::vector<size_t> inds;
			for (size_t i = cmin; i <= cmax; i++) {
				inds.push_back(i);
			}
			return inds;
		}

		friend bool operator!=( const NCPA::TridiagonalMatrix<T> &a,
				const NCPA::TridiagonalMatrix<T> &b ) {
			return a.is_unequal_( b );
		}
		friend bool operator==( const NCPA::TridiagonalMatrix<T> &a,
				const NCPA::TridiagonalMatrix<T> &b ) {
			return !(a.is_unequal_( b ));
		}

		// solves the system Ax=y for x with tridiagonal A
		std::vector<T> solve( std::vector<T> y ) const {
			T cup, pot;
			std::vector<T> diag = this->get_diagonal(),
					lower = this->get_offdiagonal(-1),
					upper = this->get_offdiagonal(1);
			size_t L = diag.size();
			if (y.size() != L) {
				throw std::runtime_error(
						"Matrix diagonal and product vector are not the same size!" );
			}
			std::vector<T> vec(L), x(L);
			int i;

			if (diag[0] == ((T)0)) {
				throw std::out_of_range( "First diagonal element is zero." );
			}
			cup = diag[0];
			vec[0] = cup;
			x[0] = y[0];
			for( i = 1; i < L; i++ ) {
				pot=lower[i-1]/cup;
				x[i]=y[i]-(x[i-1]*pot);
				cup=diag[i]-(upper[i-1]*pot);
				if ( cup == ((T)0) ) {
					throw std::out_of_range( "Matrix needs to be pivoted" );
				}
				vec[i] = cup;
			}

			x[L-1] = x[L-1] / vec[L-1];
			for ( i = L-2; i >= 0; i-- ){
				x[i] = (x[i] - (upper[i] * x[i+1]) ) / vec[i];
			}

			return x;
		}






	protected:
		bool is_unequal_( const NCPA::TridiagonalMatrix<T> &a ) const {
			if (this->rows() != a.rows()) {
				return true;
			}
			if (this->columns() != a.columns()) {
				return true;
			}
			if (this->get_diagonal() != a.get_diagonal()) {
				return true;
			}
			if (this->get_offdiagonal(-1) != a.get_offdiagonal(-1)) {
				return true;
			}
			if (this->get_offdiagonal(1) != a.get_offdiagonal(1)) {
				return true;
			}
			return false;
		}

		void check_off_diagonal_( size_t row, size_t col ) const {
			size_t mincol, maxcol;
			size_t ndiag = this->calculate_diagonal_length_();
			if (row == 0) {
				mincol = 0;
				maxcol = 1;
			} else if (row == (ndiag-1)) {
				mincol = row-1;
				maxcol = row;
			} else {
				mincol = row-1;
				maxcol = row+1;
			}
//			std::cout << "row == " << row << ", col == " << col
//				<< ", mincol == " << mincol << ", maxcol == " << maxcol << std::endl;
			if (col < mincol || col > maxcol) {
				std::ostringstream oss;
				oss << "Element [" << row << "," << col << "] is too far off-diagonal.";
				throw std::out_of_range( oss.str() );
			}
		}

		void get_column_range_( size_t rownum, size_t &cmin, size_t &cmax ) const {
			cmin = rownum > 0 ? rownum-1 : 0;
			cmax = rownum < (this->rows()-1) ? rownum+1 : this->rows()-1;
		}

		void get_row_range_( size_t colnum, size_t &rmin, size_t &rmax ) const {
			rmin = colnum > 0 ? colnum-1 : 0;
			rmax = colnum < (this->columns()-1) ? colnum+1 : this->columns()-1;
		}

	};
}

template<typename T>
void swap( NCPA::TridiagonalMatrix<T> &a, NCPA::TridiagonalMatrix<T> &b ) noexcept {
	using std::swap;
	::swap( static_cast<NCPA::SparseMatrix<T> &>(a),
			static_cast<NCPA::SparseMatrix<T> &>(b) );
}





#endif
/*
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
*/
