



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

