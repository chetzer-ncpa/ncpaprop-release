#include "Matrix.h"
#include "util.h"
#include <map>
#include <stdexcept>
#include <sstream>

namespace NCPA { template<typename T> class SparseMatrix; }

template<typename T>
void swap( NCPA::Matrix<T> &a, NCPA::Matrix<T> &b ) noexcept;

namespace NCPA {
	template<typename T>
	class SparseMatrix : public NCPA::Matrix<T> {
	public:


		// constructors
		SparseMatrix() : NCPA::Matrix() {}
		SparseMatrix( size_t nrows, size_t ncols ) : NCPA::Matrix() {
			this->initialize( nrows, ncols );
		}
		SparseMatrix( const NCPA::SparseMatrix &other ) : NCPA::Matrix() {
			this->rows_ = other.rows_;
			this->ncols_ = other.ncols_;
		}
		SparseMatrix( NCPA::SparseMatrix &&other ) noexcept : NCPA::Matrix() {
			::swap(*this,other);
		}
		virtual ~SparseMatrix() {}

		//swap
		friend void swap( NCPA::SparseMatrix &a, NCPA::SparseMatrix &b ) noexcept {
			using std::swap;
			swap( static_cast<NCPA::Matrix&>(a), static_cast<NCPA::Matrix&>(b) );
			swap( a.rows_, b.rows_ );
			swap( a.ncols_, b.ncols_ );
		}

		// assignment operator
		NCPA::SparseMatrix& operator=( NCPA::SparseMatrix other ) {
			::swap( *this, other );
			return *this;
		}

		virtual NCPA::Matrix<T> * clone() const {
			return static_cast<NCPA::Matrix<T> *>( new NCPA::SparseMatrix( *this ) );
		}

		virtual size_t rows() const {
			return rows_.size();
		}

		virtual size_t columns() const {
			return ncols_;
		}

		// initializer
		NCPA::Matrix<T> * initialize( size_t nrows, size_t ncols ) {
			rows_.clear();
			rows_ = NCPA::Vector<T>( nrows );
			ncols_ = ncols;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual T get( size_t row, size_t col ) const {
			check_dimensions_(row,col);
			auto it = rows_[row].find(col);
			if (it == rows_[row].cend()) {
				return 0;
			} else {
				return it->second;
			}
		}

		virtual NCPA::Matrix<T>* set( size_t row, size_t col, T val ) {
			check_dimensions_( row, col );
			rows_[row][col] = val;
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* get_row( size_t rownum, size_t &n, T* &values ) const {
			check_rows_( rownum );
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
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* get_row( size_t rownum, std::vector<T> &values ) const {
			check_rows_( rownum );
			values.resize( columns() );

			for (auto it = rows_[rownum].cbegin(); it != rows_[rownum].cend(); ++it) {
				values[ it->first ] = it->second;
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Vector<T> get_row( size_t rownum ) const {
			check_rows_( rownum );
			return NCPA::Vector<T>(rows_[rownum]);;
		}

		virtual NCPA::Matrix<T>* get_row( size_t rownum, NCPA::Vector<T> &values ) const {
			values = this->get_row( rownum );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* get_column( size_t colnum, size_t &n, T* &values ) const {
			check_columns_( colnum );
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
					values[ii] = it->at(rownum);
				} else {
					values[ii] = 0.0;
				}
				ii++;
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* get_column( size_t colnum, std::vector<T> &values ) const {
			check_columns_( colnum );
			values.resize( rows() );
			size_t ii = 0;
			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
				if (it->find(colnum) != it->cend()) {
					values[ii] = it->at(rownum);
				} else {
					values[ii] = 0.0;
				}
				ii++;
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Vector<T> get_column( size_t colnum ) const {
			check_columns_( colnum );
			NCPA::Vector<T> values;
			size_t ii = 0;
			for (auto it = rows_.cbegin(); it != rows_.cend(); ++it) {
				if (it->find(colnum) != it->cend()) {
					values[ii++] = it->at(rownum);
				} else {
					values[ii++] = 0.0;
				}
			}
			return values;
		}

		virtual NCPA::Matrix<T>* get_column( size_t colnum, NCPA::Vector<T> &values ) const {
			this->get_column( colnum, values );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual std::vector<T> get_diagonal() const {
			size_t n = rows() < columns() ? rows() : columns();
			std::vector<T> diag(n);
			for (size_t ii = 0; ii < n; ii++) {
				diag[n] = this->get(n,n);
			}
		}

		virtual NCPA::Matrix<T>* get_diagonal( size_t &n, T* &diag ) const {
			size_t nn = rows() < columns() ? rows() : columns();
			if (diag == nullptr) {
				n = nn;
				diag = NCPA::zeros<T>( n );
			} else {
				if (n != nn) {
					std::ostringstream oss;
					oss << "Matrix diagonal size mismatch: diagonal has "
							<< nn << " elements, but " << n << " requested.";
					throw std::out_of_range( oss.str() );
				}
			}
			for (size_t ii = 0; ii < n; ii++) {
				diag[n] = this->get(n,n);
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum, size_t nvals,
						const size_t *column_indices, const T* values ) {

			check_rows_( rownum );
			check_columns_( nvals );

			for (size_t i = 0; i < nvals; i++) {
				this->set( rownum, column_indices[i], values[i] );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const std::vector<size_t> column_indices,
						const std::vector<T> values ) {
			check_rows_( rownum );
			check_columns_( column_indices.size() );
			for (size_t i = 0; i < column_indices.size(); i++) {
				this->set( rownum, column_indices[i], values[i] );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_row( size_t rownum,
						const NCPA::Vector<T> values ) {
			check_rows_( rownum );
			check_columns_( values.size() );
			rows[ rownum ] = NCPA::Vector<T>( values );
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) {
			check_columns_( colnum );
			check_rows_( nvals );
			for (size_t i = 0; i < nvals; i++) {
				this->set( row_indices[i], colnum, values[i] );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> row_indices,
				const std::vector<T> values ) {
			check_columns_( colnum );
			check_rows_( row_indices.size() );
			for (size_t i = 0; i < row_indices.size(); i++) {
				this->set( row_indices[i], colnum, values[i] );
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T>* set_column( size_t colnum,
				const NCPA::Vector<T> values ) {
			check_columns_( colnum );
			check_rows_( values.size() );
			for (auto it = values.cbegin(); it != values.cend(); ++it) {
				this->set( it->first, colnum, it->second );
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
				NCPA::Vector<T> row = second->get_row(i);
				for (auto it = row.cbegin(); it != row.cend(); ++it) {
					this->rows_[i][it->first] += it->second;
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> *subtract( const NCPA::Matrix<T> *second ) {
			this->check_sizes_match_( second );
			for (size_t i = 0; i < second->rows(); i++) {
				NCPA::Vector<T> row = second->get_row(i);
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
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual NCPA::Matrix<T> *multiply( const NCPA::Matrix<T> *second,
				NCPA::Matrix<T> *&product ) const {
			this->check_sizes_match_transpose_( second );
			product->initialize( this->rows(), second->columns() );
			for (size_t r = 0; r < this->rows(); r++) {
				for (size_t c = 0; c < second->columns(); c++) {
					NCPA::Vector row = this->get_row(r);
					NCPA::Vector col = second->get_column(c);

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
			return product;
		}

		virtual NCPA::Matrix<T> *scale( T factor ) {
			for (auto rit = rows_.begin(); rit != rows_.end(); ++rit) {
				for (auto cit = rit->begin(); cit != rit->end(); ++cit) {
					cit->second *= factor;
				}
			}
			return static_cast<NCPA::Matrix<T> *>( this );
		}

		virtual const std::vector<size_t> &get_column_indices(size_t row) {
			if (column_indices_by_row_.size() == 0) {
				this->find_column_indices_by_row_();
			}
			return column_indices_by_row_[row];
		}

	protected:
		virtual void find_column_indices_by_row_() {
			column_indices_by_row_.clear();
			column_indices_by_row_.resize( rows_.size() );
			auto rit = rows_.cbegin();
			auto indrowit = column_indices_by_row_.begin();
			for ( ;
					rit != rows_.cend() && indrowit != column_indices_by_row_.end();
					++rit, ++indrowit) {
				indrowit->resize( rit->size() );
				for (auto cit = rit->cbegin(); cit != rit->cend(); ++cit) {
					indrowit->push_back( cit->first );
				}
			}
		}

		virtual bool is_unequal_( const NCPA::Matrix<T> &other ) const {
			if (this->rows() != other->rows()) {
				return true;
			}
			if (this->columns() != other->columns()) {
				return true;
			}
			for (size_t nrow = 0; nrow < rows(); nrow++) {
				std::vector<size_t> u_inds(this->columns());
				std::vector<size_t>::iterator ind_it =
						std::set_union(
								this->get_column_indices(nrow).cbegin(),
								this->get_column_indices(nrow).cend(),
								other.get_column_indices(nrow).cbegin(),
								other.get_column_indices(nrow).cend(),
								u_inds.begin() );
				u_inds.resize( ind_it - u_inds.begin() );
				for (ind_it = u_inds.begin();
						int_id != u_inds.end();
						++ind_it) {
					if (this->get(nrow,*ind_it) != other.get(nrow,*ind_it)) {
						return true;
					}
				}
			}
			return false;
		}


		std::vector<NCPA::Vector<T>> rows_;
		std::vector<std::vector<size_t>> column_indices_by_row_;
		size_t ncols_ = 0;

	};
}
