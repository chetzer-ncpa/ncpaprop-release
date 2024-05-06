



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
