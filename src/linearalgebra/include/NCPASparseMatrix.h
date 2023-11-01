#ifndef NCPA__LINEARALGEBRA_SPARSEMATRIX_H_INCLUDED_
#define NCPA__LINEARALGEBRA_SPARSEMATRIX_H_INCLUDED_

#include "NCPACommon.h"
#include <map>
#include <stdexcept>
#include <sstream>
#include "NCPAMatrix.h"
#include "NCPASparseVector.h"

namespace NCPA { template<typename T> class SparseMatrix; }

template<typename T>
void swap( NCPA::SparseMatrix<T> &a, NCPA::SparseMatrix<T> &b ) noexcept;
template<typename T>
bool operator==(const NCPA::SparseMatrix<T> &a, const NCPA::SparseMatrix<T> &b);
template<typename T>
bool operator!=(const NCPA::SparseMatrix<T> &a, const NCPA::SparseMatrix<T> &b);

namespace NCPA {
	template<typename T>
	class SparseMatrix : public NCPA::Matrix<T> {
	public:


		// constructors
		SparseMatrix() : NCPA::Matrix<T>() {}
		SparseMatrix( size_t nrows, size_t ncols ) : NCPA::Matrix<T>() {
			this->initialize( nrows, ncols );
		}
		SparseMatrix( const NCPA::SparseMatrix<T> &other ) : NCPA::Matrix<T>() {
			this->rows_ = other.rows_;
			this->ncols_ = other.ncols_;
		}
		SparseMatrix( NCPA::SparseMatrix<T> &&other ) noexcept : NCPA::Matrix<T>() {
			::swap(*this,other);
		}
		virtual ~SparseMatrix() {}

		//swap
		friend void ::swap<T>( NCPA::SparseMatrix<T> &a, NCPA::SparseMatrix<T> &b ) noexcept;

		// assignment operator
		NCPA::SparseMatrix<T>& operator=( NCPA::SparseMatrix<T> other ) {
			swap( *this, other );
			return *this;
		}

		virtual NCPA::Matrix<T> * clone() const {
			return static_cast<NCPA::Matrix<T> *>( new NCPA::SparseMatrix<T>( *this ) );
		}

		virtual size_t rows() const {
			return rows_.size();
		}

		virtual size_t columns() const {
			return ncols_;
		}

		// initializer
		virtual NCPA::Matrix<T> * initialize( size_t nrows, size_t ncols ) {
			rows_.clear();
			rows_ = std::vector<NCPA::SparseVector<T>>( nrows );
			ncols_ = ncols;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual T get( size_t row, size_t col ) const {
			this->check_dimensions(row,col);
			auto it = rows_[row].find(col);
			if (it == rows_[row].cend()) {
				return 0;
			} else {
				return it->second;
			}
		}

		virtual NCPA::Matrix<T>* set( size_t row, size_t col, T val ) {
			this->check_dimensions( row, col );
			rows_[row][col] = val;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual void get_row( size_t rownum, size_t &n, T* &values ) const {
			this->check_rows_( rownum );
			if (values == nullptr) {
				n = columns();
				values = NCPA::zeros<T>( n );
			} else if (n != columns()) {
				std::ostringstream oss;
				oss << "Matrix row size mismatch: Matrix has " << rows()
						<< " rows, but " << n << " values requested.";
				throw std::out_of_range( oss.str() );
			}
			for (auto it = rows_[rownum].cbegin(); it != rows_[rownum].cend(); ++it) {
				values[ it->first ] = it->second;
			}
		}

		virtual void get_row( size_t rownum, std::vector<T> &values ) const {
			this->check_rows_( rownum );
			values.resize( columns() );

			for (auto it = rows_[rownum].cbegin(); it != rows_[rownum].cend(); ++it) {
				values[ it->first ] = it->second;
			}
		}

		virtual NCPA::SparseVector<T> get_row( size_t rownum ) const {
			this->check_rows_( rownum );
			return NCPA::SparseVector<T>(rows_[rownum]);
		}

		virtual void get_row( size_t rownum, NCPA::SparseVector<T> &values ) const {
			values = this->get_row( rownum );
		}

		virtual void get_column( size_t colnum, size_t &n, T* &values ) const {
			this->check_columns_( colnum );
			if (values == nullptr) {
				n = rows();
				values = NCPA::zeros<T>( n );
			} else if (n != rows()) {
				std::ostringstream oss;
				oss << "Matrix column size mismatch: Matrix has " << columns()
						<< " columns, but " << n << " values requested.";
				throw std::out_of_range( oss.str() );
			}
			size_t ii = 0;
			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
				if (it->find(colnum) != it->cend()) {
					values[ii] = it->at(colnum);
				} else {
					values[ii] = 0.0;
				}
				ii++;
			}
		}

		virtual void get_column( size_t colnum, std::vector<T> &values ) const {
			this->check_columns_( colnum );
			values.resize( rows() );
			size_t ii = 0;
			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
				if (it->find(colnum) != it->cend()) {
					values[ii] = it->at(colnum);
				} else {
					values[ii] = 0.0;
				}
				ii++;
			}
		}

		virtual NCPA::SparseVector<T> get_column( size_t colnum ) const {
			this->check_columns_( colnum );
			NCPA::SparseVector<T> values;
			size_t ii = 0;
			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
				if (it->find(colnum) != it->cend()) {
					values[ii] = it->at(colnum);
//				} else {
//					values[ii++] = 0.0;
				}
				ii++;
			}
			return values;
		}

		virtual void get_column( size_t colnum, NCPA::SparseVector<T> &values ) const {
			values = this->get_column( colnum );
		}


		virtual NCPA::Matrix<T>* set_row( size_t rownum, size_t nvals,
						const size_t * column_indices, const T* values ) {

			this->check_rows_( rownum );
			this->check_columns_( nvals );
			rows_[rownum].clear();
			for (size_t i = 0; i < nvals; i++) {
				rows_[rownum][column_indices[i]] = values[i];
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<T> &values ) {
			this->check_rows_( rownum );
			this->check_columns_( values.size() );
			rows_[rownum].clear();
			for (size_t i = 0; i < values.size(); i++) {
				if (values[i] != (T)0.0) {
					rows_[rownum][i] = values[i];
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<size_t> &column_indices,
						const std::vector<T> &values ) {
			this->check_rows_( rownum );
			this->check_columns_( column_indices.size() );
			rows_[rownum].clear();
			for (size_t i = 0; i < values.size(); i++) {
				rows_[rownum][column_indices[i]] = values[i];
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const NCPA::SparseVector<T> &values ) {
			this->check_rows_( rownum );
			this->check_columns_( values.size() );
			rows_[ rownum ] = NCPA::SparseVector<T>( values );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) {
			this->check_columns_( colnum );
			this->check_rows_( nvals );
			for (size_t i = 0; i < rows(); i++) {
				auto cit = rows_[i].find( colnum );
				if (cit != rows_[i].end()) {
					rows_[i].erase(cit);
				}
			}
			for (size_t i = 0; i < nvals; i++) {
//				std::cout << "Setting rows_[" << row_indices[i] << "][" << colnum << "] = "
//						<< values[i] << std::endl;
				rows_[row_indices[i]][colnum] = values[i];
//				this->set( row_indices[i], colnum, values[i] );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<T> &values ) {
			this->check_columns_( colnum );
			this->check_rows_( values.size() );
			for (size_t i = 0; i < rows(); i++) {
				if (values[i] != (T)0) {
					rows_[i][colnum] = values[i];
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> &row_indices,
				const std::vector<T> &values ) {
			this->check_columns_( colnum );
			this->check_rows_( row_indices.size() );
			for (size_t i = 0; i < rows(); i++) {
				auto cit = rows_[i].find( colnum );
				if (cit != rows_[i].end()) {
					rows_[i].erase(cit);
				}
			}
			for (size_t i = 0; i < values.size(); i++) {
				rows_[row_indices[i]][colnum] = values[i];
//				this->set( row_indices[i], colnum, values[i] );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const NCPA::SparseVector<T> &values ) {
			this->check_columns_( colnum );
			this->check_rows_( values.size() );
			for (size_t i = 0; i < rows(); i++) {
				auto cit = rows_[i].find( colnum );
				if (cit != rows_[i].end()) {
					rows_[i].erase(cit);
				}
			}
			for (auto it = values.cbegin(); it != values.cend(); ++it) {
				rows_[it->first][colnum] = it->second;
//				this->set( it->first, colnum, it->second );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}




		virtual NCPA::Matrix<T>* ready() {
			// possibly do some caching in here
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual bool is_ready() const {
			return (ncols_ > 0) && (rows_.size() > 0);
		}

		virtual NCPA::Matrix<T> *add( const NCPA::Matrix<T> *second ) {
			this->check_sizes_match_( second );
			for (size_t i = 0; i < second->rows(); i++) {
				NCPA::SparseVector<T> row;
				second->get_row(i,row);
				for (auto it = row.cbegin(); it != row.cend(); ++it) {
					this->rows_[i][it->first] += it->second;
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> *subtract( const NCPA::Matrix<T> *second ) {
			this->check_sizes_match_( second );
			for (size_t i = 0; i < second->rows(); i++) {
				NCPA::SparseVector<T> row;
				second->get_row(i, row);
				for (auto it = row.cbegin(); it != row.cend(); ++it) {
					this->rows_[i][it->first] -= it->second;
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> *transpose() {
			NCPA::SparseMatrix<T> other( this->columns(), this->rows() );
			size_t counter = 0;
			for (auto rit = rows_.begin(); rit != rows_.end(); ++rit) {
				other.set_column(counter++, *rit);
			}
			swap( *this, other );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		static void multiply(
				const NCPA::Matrix<T> *first,
				const NCPA::SparseMatrix<T> *second,
				NCPA::Matrix<T> *&product ) {
			NCPA::Matrix<T>::multiply(
					first,
					static_cast<NCPA::Matrix<T>*>( second ),
					product );
		}

		static void multiply(
				const NCPA::SparseMatrix<T> *first,
				const NCPA::Matrix<T> *second,
				NCPA::Matrix<T> *&product ) {
			NCPA::Matrix<T>::multiply(
					static_cast<NCPA::Matrix<T>*>( first ),
					second,
					product );
		}

		static void multiply(
				const NCPA::SparseMatrix<T> *first,
				const NCPA::SparseMatrix<T> *second,
				NCPA::Matrix<T> *&product ) {
			if (first->columns() != second->rows()) {
				std::ostringstream oss;
				oss << "Size mismatch in matrix product: Columns in first matrix ("
						<< first->columns() << ") must equal rows in second matrix ("
						<< second->rows() << ")";
				throw std::out_of_range( oss.str() );
			}
			if (product == nullptr) {
				product = static_cast<NCPA::Matrix<T> *>(
						new NCPA::SparseMatrix<T>( first->rows(), second->columns() )
					);
			} else {
				product->initialize( first->rows(), second->columns() );
			}

			for (size_t r = 0; r < first->rows(); r++) {
				for (size_t c = 0; c < second->columns(); c++) {
					NCPA::SparseVector<T> row, col;
					first->get_row(r, row);
					second->get_column(c, col);

					auto rit = row.cbegin();
					auto cit = col.cbegin();
					T prod = 0;
					while (rit != row.cend() || cit != col.cend()) {
						// product will only be nonzero if the indices match
						if (rit->first == cit->first) {
							prod += rit->second * cit->second;
						} else if (rit->first > cit->first) {
							++cit;
						} else {
							++rit;
						}
					}
					product->set( r, c, prod );
				}
			}
		}

		virtual NCPA::Matrix<T> *scale( T factor ) {
			for (auto rit = rows_.begin(); rit != rows_.end(); ++rit) {
				for (auto cit = rit->begin(); cit != rit->end(); ++cit) {
					cit->second *= factor;
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual std::vector<size_t> get_column_indices(size_t row) const {
			assert_ready_();
			std::vector<size_t> inds;
			for (auto it = rows_[row].cbegin(); it != rows_[row].cend(); ++it) {
				inds.push_back( it->first );
			}
			std::sort( inds.begin(), inds.end() );
			return inds;
		}

		// operators
		friend bool operator!=( const NCPA::SparseMatrix<T> &a, const NCPA::SparseMatrix<T> &b ) {
			return a.is_unequal_( b );
		}
		friend bool operator==( const NCPA::SparseMatrix<T> &a, const NCPA::SparseMatrix<T> &b ) {
			return !(a.is_unequal_( b ));
		}


	protected:
		virtual void assert_ready_() const {
			if (!this->is_ready()) {
				throw std::logic_error( "Matrix not ready!" );
			}
		}

		virtual bool is_unequal_( const NCPA::Matrix<T> &other ) const {
			if (this->rows() != other.rows()) {
//				std::cout << "Row mismatch" << std::endl;
				return true;
			}
			if (this->columns() != other.columns()) {
//				std::cout << "Row mismatch" << std::endl;
				return true;
			}
//			std::cout << "Checking values" << std::endl;
			for (size_t nrow = 0; nrow < rows(); nrow++) {
				std::vector<size_t> u_inds(2*this->columns()),
						thiscols = this->get_column_indices(nrow),
						othercols = other.get_column_indices(nrow);
				std::vector<size_t>::iterator ind_it =
						std::set_union(
								thiscols.cbegin(),
								thiscols.cend(),
								othercols.cbegin(),
								othercols.cend(),
								u_inds.begin() );
				u_inds.resize( ind_it - u_inds.begin() );
				for (ind_it = u_inds.begin();
						ind_it != u_inds.end();
						++ind_it) {
					if (this->get(nrow,*ind_it) != other.get(nrow,*ind_it)) {
//						std::cout << "Mismatch in row " << nrow << std::endl;
						return true;
					}
				}
			}
			return false;
		}


		std::vector<NCPA::SparseVector<T>> rows_;
		size_t ncols_ = 0;

	};
}

template<typename T>
void swap( NCPA::SparseMatrix<T> &a, NCPA::SparseMatrix<T> &b ) noexcept {
	using std::swap;
	::swap( static_cast<NCPA::Matrix<T> &>(a), static_cast<NCPA::Matrix<T> &>(b) );
	swap( a.rows_, b.rows_ );
	swap( a.ncols_, b.ncols_ );
}


#endif
