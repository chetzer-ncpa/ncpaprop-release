#ifndef NCPA_LINEARALGEBRA__NCPALINEARSYSTEMSOLVER_H_INCLUDED
#define NCPA_LINEARALGEBRA__NCPALINEARSYSTEMSOLVER_H_INCLUDED

#include <stdexcept>
#include "NCPAMatrix.h"
#include "NCPAVector.h"
#include "NCPALUDecomposition.h"

namespace NCPA {
	template<typename T> class LinearSystemSolver;

	enum class solver_method_t : size_t {
		NONE = 0,

		CHOOSE_BEST,

		// Solve using LU decomposition
		LU_DECOMPOSITION_PIVOT,
		LU_DECOMPOSITION_NO_PIVOT,

		// Solve using efficient tridiagonal algorithm
		TRIDIAGONAL,
	};
}

template<typename T>
void swap( NCPA::LinearSystemSolver<T> &a,
		NCPA::LinearSystemSolver<T> &b ) noexcept;

namespace NCPA {



	template<typename T> class LinearSystemSolver {

	public:
		LinearSystemSolver() = default;
		LinearSystemSolver( const NCPA::LinearSystemSolver<T> &other ) {
			if (other.mat_ != nullptr) {
				this->mat_ = other.mat_->clone();
			}
			if (other.decomp_ != nullptr) {
				this->decomp_ = other.decomp_->clone();
			}
			this->method_ = other.method_;
		}
		LinearSystemSolver( NCPA::LinearSystemSolver<T> &&other ) noexcept {
			::swap( *this, other );
		}

		virtual ~LinearSystemSolver() {
			this->clear();
		}

		friend void swap<T>( NCPA::LinearSystemSolver<T> &a,
				NCPA::LinearSystemSolver<T> &b ) noexcept;

		// resets internals
		virtual NCPA::LinearSystemSolver<T> * set_system_matrix(
				NCPA::Matrix<T> * m ) = 0;
		virtual NCPA::LinearSystemSolver<T> * compute_lu_decomposition(
				bool pivot = true ) = 0;
		virtual NCPA::solver_method_t choose_best_method() const = 0;


		virtual NCPA::LinearSystemSolver<T> * set_method(
				NCPA::solver_method_t method ) {
			if (method == NCPA::solver_method_t::CHOOSE_BEST) {
				method_ = this->choose_best_method();
			} else {
				method_ = method;
			}
			return this;
		}

		virtual NCPA::LinearSystemSolver<T> * clear() {
			if (decomp_ != nullptr) {
				delete decomp_;
				decomp_ = nullptr;
			}
			if (mat_ != nullptr) {
				delete mat_;
				mat_ = nullptr;
			}
			return this;
		}

		virtual NCPA::LinearSystemSolver<T> * solve(
				const NCPA::Vector<T> *b_in,
				NCPA::Vector<T> *x_out ) {
			switch (method_) {
			case solver_method_t::NONE:
				throw std::invalid_argument( "No solution method specified!" );
				break;
			case solver_method_t::LU_DECOMPOSITION_PIVOT:
				return this->solve_lu( b_in, x_out, true );
				break;
			case solver_method_t::LU_DECOMPOSITION_NO_PIVOT:
				return this->solve_lu( b_in, x_out, false );
				break;
			case solver_method_t::TRIDIAGONAL:
				return this->solve_tridiagonal( b_in, x_out );
				break;
			default:
				throw std::out_of_range(
						"Unknown or unimplemented solver type requested.");
			}
			return nullptr;
		}

		virtual NCPA::LinearSystemSolver<T> * solve_lu(
				const NCPA::Vector<T> *b,
				NCPA::Vector<T> *x,
				bool pivot = true ) {

			if (this->decomp_ == nullptr) {
				this->compute_lu_decomposition( pivot );
			}

			// make sure vector size is correct
			this->decomp_->lower()->assert_rows_match( b.size() );

			size_t N = lower->rows();
			int i=0, j=0, iN = (int)N;

			// temporary vectors
			std::vector<T> y( N, ((T)0) ), Pb( N, ((T)0) );

			// forward substitution
			for (i = 0; i < iN; i++) {

				// @todo is this loop necessary if no pivoting?
				for (j = 0; j < iN; j++) {
					Pb[i] +=
							(this->decomp_->permutation()->get(i,j))
							* b->get(j);
				}
				for (j = 0; j < i; j++) {
					Pb[i] -= this->decomp_->lower()->get(i,j) * y[j];
				}
				y[i] = Pb[i] / this->decomp_->lower()->get(i,i);
			}

			for (i = iN-1; i >= 0; i--) {
				for (j = i+1; j < iN; j++) {
					y[i] -= this->decomp_->upper()->get(i,j)
							* x->get(j);
				}
				x->set( i, y[i] / upper->get(i,i) )
			}

			return this;
		}

		virtual NCPA::LinearSystemSolver<T> * solve_tridiagonal(
				const NCPA::Vector<T> *b_in,
				NCPA::Vector<T> *x_out ) {

			if (!this->mat_->is_square()) {
				throw std::domain_error(
						"solve_tridiagonal() not implemented for non-square matrices" );
			}
			T cup, pot;
			std::vector<T> diag = this->get_diagonal(),
					lower = this->get_offdiagonal(-1),
					upper = this->get_offdiagonal(1),
					b = b_in->as_std();

			size_t L = diag.size();
			if (b.size() != L) {
				throw std::invalid_argument(
						"Matrix diagonal and RHS vector are not the same size!" );
			}
			if (b->size() != L) {
				throw std::invalid_argument(
						"Matrix diagonal and solution vector are not the same size!" );
			}
			std::vector<T> vec(L), x(L);
			int i;

			if (diag[0] == ((T)0)) {
				throw std::out_of_range( "First diagonal element is zero." );
			}
			cup = diag[0];
			vec[0] = cup;
			x[0] = b[0];
			for( i = 1; i < L; i++ ) {
				pot=lower[i-1]/cup;
				x[i]=b[i]-(x[i-1]*pot);
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
			x_out->from_std( x );

			return this;
		}

	protected:
		NCPA::LUDecomposition<T> *decomp_ = nullptr;
		NCPA::Matrix<T> *mat_ = nullptr;
		NCPA::solver_method_t method_ = NCPA::solver_method_t::NONE;
	};
}

template<typename T>
void swap( NCPA::LinearSystemSolver<T> &a,
		NCPA::LinearSystemSolver<T> &b ) noexcept {
	using std::swap;
	swap( a.decomp_, b.decomp_ );
	swap( a.mat_, b.mat_ );
	swap( a.method_, b.method_ );
}

#endif
