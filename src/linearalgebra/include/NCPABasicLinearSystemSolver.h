#ifndef NCPA_LINEARALGEBRA__NCPABASICLINEARSYSTEMSOLVER_H_INCLUDED
#define NCPA_LINEARALGEBRA__NCPABASICLINEARSYSTEMSOLVER_H_INCLUDED

#include <stdexcept>
#include "NCPALinearSystemSolver.h"
#include "NCPAMatrix.h"
#include "NCPABasicMatrix.h"
#include "NCPAVector.h"
#include "NCPABasicVector.h"
#include "NCPABasicLUDecomposition.h"


namespace NCPA {
	template<typename T> class BasicLinearSystemSolver
			: public NCPA::LinearSystemSolver<T> {

	public:
		BasicLinearSystemSolver() = default;
		virtual ~BasicLinearSystemSolver() {}

		// resets internals
		virtual NCPA::LinearSystemSolver<T> * set_system_matrix( NCPA::Matrix<T> * m ) {
			this->clear();
			this->mat_ = new NCPA::BasicMatrix<T>();
			this->mat_->copy( m );
			return static_cast<NCPA::LinearSystemSolver<T>*>( this );
		}

		virtual NCPA::LinearSystemSolver<T> * compute_lu_decomposition(
				bool pivot = true ) {
			if (this->decomp_ != nullptr) {
				delete this->decomp_;
			}
			if (this->mat_ == nullptr) {
				throw std::invalid_argument( "No system matrix set!" );
			}
			this->decomp_ = new NCPA::BasicLUDecomposition<T>();
			this->decomp_->initialize( this->mat_ );
			this->decomp_->compute( this->mat_, pivot );
			return static_cast<NCPA::LinearSystemSolver<T> *>( this );
		}

		virtual NCPA::solver_method_t choose_best_method() const {
			if (this->mat_ == nullptr) {
				throw std::invalid_argument( "Cannot determine best solution method for null Matrix.");
			}
			if (this->mat_->is_tridiagonal()) {
				return NCPA::solver_method_t::TRIDIAGONAL;
			} else {
				return NCPA::solver_method_t::LU_DECOMPOSITION_PIVOT;
			}
		}
	};



#endif
