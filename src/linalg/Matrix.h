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
#include <unordered_map>
#include "util.h"



namespace NCPA { template<typename T> class Matrix; }

template<typename T>
void swap( NCPA::Matrix<T> &a, NCPA::Matrix<T> &b ) noexcept;

template<typename T>
bool operator==(const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &a);
template<typename T>
bool operator!=(const NCPA::Matrix<T> &a, const NCPA::Matrix<T> &a);

namespace NCPA {

	// Vector template, used for rows and columns as well
	template<typename T>
	using Vector = std::unordered_map<size_t,T>;

	template<typename T>
	class Matrix {

	public:
		Matrix() {}
		Matrix( const Matrix &other ) {}
		Matrix( Matrix &&other ) noexcept {}
		virtual ~Matrix() {}
		friend void ::swap( Matrix<T> &a, Matrix<T> &b ) noexcept {}
		virtual Matrix<T>* clone() const = 0;

		// element access
		virtual size_t rows() const = 0;
		virtual size_t columns() const = 0;
		virtual T get( size_t d1, size_t d2 ) const = 0;
		virtual Matrix<T>* set( size_t d1, size_t d2, T val ) = 0;


		// row and column access
		virtual Matrix<T>* get_row( size_t rownum, size_t &n, T* &values ) = 0;
		virtual Matrix<T>* get_row( size_t rownum, std::vector<T> &values ) = 0;
		virtual Matrix<T>* get_row( size_t rownum, NCPA::Vector<T> &row ) = 0;
		virtual NCPA::Vector<T> get_row( size_t rownum ) = 0;
		virtual Matrix<T>* get_column( size_t colnum, size_t &n, T* &values ) = 0;
		virtual Matrix<T>* get_column( size_t colnum, std::vector<T> &values ) = 0;
		virtual Matrix<T>* get_column( size_t colnum, NCPA::Vector<T> &row ) = 0;
		virtual NCPA::Vector<T> get_column( size_t colnum ) = 0;
		virtual std::vector<T> get_diagonal() const = 0;
		virtual Matrix<T>* get_diagonal( size_t &n, T* &vals ) const = 0;

		virtual Matrix<T>* set_row( size_t rownum, size_t nvals,
				const size_t *column_indices, const T* values ) = 0;
		virtual Matrix<T>* set_row( size_t rownum,
				const std::vector<size_t> column_indices,
				const std::vector<T> values ) = 0;
		virtual Matrix<T>* set_row( size_t rownum,
				const NCPA::Vector<T> values ) = 0;

		virtual Matrix<T>* set_column( size_t colnum, size_t nvals,
				const size_t *row_indices, const T* values ) = 0;
		virtual Matrix<T>* set_column( size_t colnum,
				const std::vector<size_t> row_indices,
				const std::vector<T> values ) = 0;
		virtual Matrix<T>* set_column( size_t colnum,
				const NCPA::Vector<T> values ) = 0;

		NCPA::Vector<T>& operator[]( int rownum ) {
			NCPA::Vector<T> row;
			this->get_row( rownum, row );
			return row;
		}

		// status
		virtual Matrix<T>* ready() = 0;
		virtual bool is_ready() const = 0;

		// operations
		virtual Matrix<T> *add( const Matrix<T> *second ) = 0;
		virtual Matrix<T> *subtract( const Matrix<T> *second ) = 0;
		virtual Matrix<T> *multiply( const Matrix<T> *second, Matrix<T> *&product ) = 0;
		virtual Matrix<T> *transpose() = 0;
		virtual Matrix<T> *scale( T factor ) = 0;

		// operators
		friend bool operator!=( const Matrix<T> other ) {
			return this->is_unequal_( other );
		}
		friend bool operator==( const Matrix<T> other ) {
			return !(this->is_unequal_( other ));
		}


	protected:
		virtual bool is_unequal_( const Matrix<T> &other ) const = 0;
		virtual void check_dimensions_( size_t d1, size_t d2 ) const {
			this->check_rows_(d1);
			this->check_columns_(d2);
		}
		virtual void check_rows_( size_t d1 ) const {
			if (d1 >= this->rows()) {
				std::ostringstream oss;
				oss << "Requested row " << d1
						<< "is out of range for Matrix with "
						<< rows() << " rows.";
				throw std::out_of_range( oss.str() );
			}
		}
		virtual void check_columns_( size_t d1 ) const {
			if (d1 >= this->columns()) {
				std::ostringstream oss;
				oss << "Requested column " << d1
						<< "is out of range for Matrix with "
						<< columns() << " columns.";
				throw std::out_of_range( oss.str() );
			}
		}
		virtual void check_sizes_match_( Matrix<T> *other ) {
			if (this->rows() != other->rows() || this->columns() != other->columns()) {
				std::ostringstream oss;
				oss << "Matrix size mismatch between RHS matrix of size "
						<< other->rows() << "x" << other->columns()
						<< " and LHS matrix of size " << this->rows() << "x"
						<< this->columns() << std::endl;
				throw std::out_of_range( oss.str() );
			}
		}
		virtual void check_sizes_match_transpose_( Matrix<T> *other ) {
			if (this->rows() != other->columns() || this->columns() != other->rows()) {
				std::ostringstream oss;
				oss << "Matrix transpose size mismatch between RHS matrix of size "
						<< other->columns() << "x" << other->rows()
						<< " and LHS matrix of size " << this->rows() << "x"
						<< this->columns() << std::endl;
				throw std::out_of_range( oss.str() );
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
