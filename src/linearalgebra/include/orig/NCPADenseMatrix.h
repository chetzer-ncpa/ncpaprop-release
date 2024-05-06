


	template<typename T>
	class DenseMatrix : public Matrix<T> {

	protected:
		T **contents_;
		size_t rows_, cols_;
		bool ready_;

		virtual bool check_dimensions( size_t d1, size_t d2 ) const {
			return ( d1 < rows_ ) && ( d2 < cols_ );
		}

	public:
		DenseMatrix( size_t d1, size_t d2 ) {
			rows_ = (size_t)d1;
			cols_ = (size_t)d2;
			contents_ = NCPA::allocate_matrix<T>( rows_, cols_ );
			ready();
		}

		DenseMatrix( const NCPA::DenseMatrix<T> &in ) {
			rows_ = in.rows();
			cols_ = in.columns();
			contents_ = NCPA::allocate_matrix<T>( rows_, cols_ );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					contents_[i][j] = in.get( i, j );
				}
			}
			ready();
		}

		virtual ~DenseMatrix() {
			NCPA::free_matrix<T>( contents_, rows_, cols_ );
		}

		virtual T get( size_t d1, size_t d2 ) const {
			assert(ready_);
			check_dimensions( d1, d2 );
			return contents_[ d1 ][ d2 ];
		}

		virtual T at( size_t d1, size_t d2 ) const {
			return this->get( d1, d2 );
		}

		virtual void set( size_t d1, size_t d2, T val ) {
			check_dimensions( d1, d2 );
			contents_[ d1 ][ d2 ] = val;
		}

		// nothing special needs to be done to make ready
		virtual void ready() {
			ready_ = true;
		}

		virtual bool is_ready() const {
			return ready_;
		}

		virtual size_t rows() const {
			return rows_;
		}

		virtual size_t columns() const {
			return cols_;
		}

		virtual Matrix<T> *add( const Matrix<T> *term ) const {
			assert( this->rows() == term->rows() );
			assert( this->columns() == term->columns() );
			assert( ready_ );
			assert( term->is_ready() );

			Matrix<T> *sum = new DenseMatrix<T>( rows_, cols_ );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					sum->set( i, j, contents_[ i ][ j ] + term->get( i, j ) );
				}
			}

			sum->ready();
			return sum;
		}

		virtual Matrix<T> *multiply( const Matrix<T> *second ) const {
			assert( this->columns() == second->rows() );
			assert( ready_ );
			assert( second->is_ready() );

			Matrix<T> *product = new DenseMatrix<T>( rows_, second->columns() );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < second->columns(); j++) {
					T temp = 0;
					for (size_t k = 0; k < cols_; k++) {
						temp += contents_[i][k]
							* second->get(k,j);
					}
					product->set( i, j, temp );
				}
			}

			product->ready();
			return product;
		}

		virtual Matrix<T> *transpose() const {
			assert( ready_ );

			Matrix<T> *trans = new DenseMatrix<T>( cols_, rows_ );
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					trans->set( j, i, contents_[ i ][ j ] );
				}
			}

			trans->ready();
			return trans;
		}

		virtual void scale( T factor ) {
			for (size_t i = 0; i < rows_; i++) {
				for (size_t j = 0; j < cols_; j++) {
					contents_[ i ][ j ] *= factor;
				}
			}
		}

		virtual void get_row( size_t row, size_t &n,
				size_t *&indices, T *&values ) const {
			assert( row < rows_ );
			n = cols_;
			indices = NCPA::index_vector<size_t>( n );
			values  = NCPA::zeros<T>( n );
			std::memcpy( values, contents_[row], n*sizeof(T) );
		}

		virtual void get_row( size_t row, size_t &n,
				T *&values ) const {
			assert( row < rows_ );
			n = cols_;
			values  = NCPA::zeros<T>( n );
			std::memcpy( values, contents_[row], n*sizeof(T) );
		}

		virtual void get_row( size_t row,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const {
			assert( row < rows_ );
			indices.reserve( cols_ );
			values.reserve( cols_ );
			for (size_t i = 0; i < cols_; i++) {
				indices.push_back( i );
				values.push_back( contents_[row][i] );
			}
		}

		virtual void get_row( size_t row,
				std::vector<T> &values ) const {
			assert( row < rows_ );
			values.reserve( cols_ );
			for (size_t i = 0; i < cols_; i++) {
				values.push_back( contents_[row][i] );
			}
		}

		virtual void get_row_indices( size_t row, size_t &n,
				size_t *&indices ) const {
			assert( row < rows_ );
			n = cols_;
			indices = NCPA::index_vector<size_t>( n );
		}

		virtual void get_row_indices( size_t row,
				std::vector<size_t> &indices ) const {
			assert( row < rows_ );
			indices.reserve( cols_ );
			for (size_t i = 0; i < cols_; i++) {
				indices.push_back( i );
			}
		}


		virtual void get_column( size_t col, size_t &n,
				size_t *&indices, T *&values ) const {
			assert( col < cols_ );
			n = rows_;
			indices = NCPA::index_vector<size_t>( n );
			values  = NCPA::zeros<T>( n );
			for (size_t i = 0; i < n; i++) {
				values[ i ] = contents_[ i ][ col ];
			}
		}

		virtual void get_column( size_t col, size_t &n,
				T *&values ) const {
			assert( col < cols_ );
			n = rows_;
			values  = NCPA::zeros<T>( n );
			for (size_t i = 0; i < n; i++) {
				values[ i ] = contents_[ i ][ col ];
			}
		}

		virtual void get_column( size_t col,
				std::vector<size_t> &indices,
				std::vector<T> &values ) const {
			assert( col < cols_ );
			indices.reserve( rows_ );
			values.reserve( rows_ );
			for (size_t i = 0; i < rows_; i++) {
				indices.push_back( i );
				values.push_back( contents_[i][col] );
			}
		}

		virtual void get_column( size_t col,
				std::vector<T> &values ) const {
			assert( col < cols_ );
			values.reserve( rows_ );
			for (size_t i = 0; i < rows_; i++) {
				values.push_back( contents_[i][col] );
			}
		}

		virtual void get_column_indices( size_t col, size_t &n,
				size_t *&indices ) const {
			assert( col < cols_ );
			n = rows_;
			indices = NCPA::index_vector<size_t>( n );
		}

		virtual void get_column_indices( size_t col,
				std::vector<size_t> &indices ) const {
			assert( col < cols_ );
			indices.reserve( rows_ );
			for (size_t i = 0; i < rows_; i++) {
				indices.push_back( i );
			}
		}

		virtual size_t get_row_capacity() const {
			return cols_;
		}


	};
