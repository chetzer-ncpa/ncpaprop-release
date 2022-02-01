#ifndef NCPAPROP_UTIL_H_INCLUDED
#define NCPAPROP_UTIL_H_INCLUDED

#include <string>
#include <istream>
#include <iostream>
#include <vector>
#include <complex>
#include <cfloat>
#include <cstring>
#include <cassert>

#ifndef PI
#define PI 3.141592653589793
#endif

namespace NCPA {

	std::string timeAsString( double d );
	bool fexists( const char *filename );
	std::string deblank( const std::string& orig );

	double mean( double*, size_t );
	size_t nextpow2( size_t v );

	std::istream &safe_getline( std::istream &is, std::string &s );

	// perform a simple linear interpolation
	double linearInterp( double x1, double y1, double x2, double y2, double x );
	void toUpperCase( char *in );
	void toLowerCase( char *in );
	std::string toUpperCase( const std::string in );
	std::string toLowerCase( const std::string in );

	// Split a string into more strings by tokenizing
	std::vector< std::string > split( std::string input, std::string delimiters = " \t\n" );
	std::string deblank( const std::string& str, const std::string& whitespace );

	// random numbers
	std::vector<double> random_numbers( size_t N, double scale = 1.0 );

	bool checkAzimuthLimits( double toCheck, double target, double tolerance );
	// double normalizeAzimuth( double in );
	template<typename T> T normalizeAzimuth( T in ) {
		double out = (double)in;
		while (out < 0.0) {
			out += 360.0;
		}
		while (out >= 360.0) {
			out -= 360.0;
		}
		return (T)out;
	}

	template<typename T> T normalizeAzimuthRadians( T in ) {
		double out = (double)in;
		while (out < 0.0) {
			out += 2.0 * PI;
		}
		while (out >= 2.0*PI) {
			out -= 2.0 * PI;
		}
		return (T)out;
	}
	
	// Utility functions
	int **imatrix(long nr, long nc);
	int free_imatrix(int** v, long nr, long nc);
	double **dmatrix(long nr, long nc);
	int free_dmatrix(double** v, long nr, long nc);
	std::complex<double> **cmatrix(long nr, long nc);
	int free_cmatrix(std::complex<double> **v, long nr, long nc);
	std::complex<double> ***c3Darray(size_t xlen, size_t ylen, size_t zlen);
	void free_c3Darray(std::complex<double> ***data, size_t xlen, size_t ylen);

	// coordinate system transformations
	void pol2cart( double r, double theta_rad, double &x, double &y );
	void cart2pol( double x, double y, double &r, double &theta_rad );

	int count_rows_arbcol( const std::string& filename );
	void read_matrix_from_file( const std::string &filename, double **&contents,
		size_t &ncols, size_t &nrows, std::string delimiters = " \t\n" );
	void read_text_columns_from_file( const std::string &filename,
		std::vector< std::vector< std::string > > &contents,
		std::string delimiters = " \t\n" );

	// function templates.  Need to be defined here as well as declared so compiler can
	// find and generate them properly.
	template<typename T> T min( T a, T b ) { return a < b ? a : b; }
	template<typename T> T max( T a, T b ) { return a > b ? a : b; }

	template<typename T> T max( const T *vals, size_t size ) {
		T maxval = vals[ 0 ];
		for (size_t i = 1; i < size; i++) {
			maxval = NCPA::max( vals[i], maxval );
		}
		return maxval;
	}

	template<typename T> T min( const T *vals, size_t size ) {
		T minval = vals[ 0 ];
		for (size_t i = 1; i < size; i++) {
			minval = NCPA::min( vals[ i ], minval );
		}
		return minval;
	}

	template<typename T> T max( std::vector< T > vals ) {
		T maxval = vals.front();
		for (typename std::vector<T>::const_iterator cit = vals.cbegin();
				cit != vals.cend(); ++cit) {
			maxval = NCPA::max( *cit, maxval );
		}
		return maxval;
	}

	template<typename T> T min( std::vector< T > vals ) {
		T minval = vals.front();
		for (typename std::vector<T>::const_iterator cit = vals.cbegin();
				cit != vals.cend(); ++cit) {
			minval = NCPA::min( *cit, minval );
		}
		return minval;
	}

	template<typename T> T deg2rad( T deg_in ) { return (T)(((double)deg_in) * PI / 180.0); }
	template<typename T> T rad2deg( T deg_in ) { return (T)(((double)deg_in) * 180.0 / PI); }

	template<typename T> void move_origin( size_t npts,
			const T *old_x, const T *old_y, T x_new_origin, T y_new_origin,
			T *new_x, T *new_y ) {
		for (size_t i = 0; i < npts; i++) {
			new_x[ i ] = old_x[ i ] - x_new_origin;
			new_y[ i ] = old_y[ i ] - y_new_origin;
		}
	}

	template<typename T> void move_origin( T old_x, T old_y,
			T x_new_origin, T y_new_origin, T &new_x, T &new_y ) {
		new_x = old_x - x_new_origin;
		new_y = old_y - y_new_origin;
	}

	template<typename T> T *zeros( size_t n ) {
		T *out = new T[ n ];
		std::memset( out, 0, n*sizeof(T) );
		return out;
	}

	template<typename T> T *single_valued_vector( size_t n, T val ) {
		T *out = new T[ n ];
		for (size_t i = 0; i < n; i++) {
			out[ i ] = val;
		}
		return out;
	}

	template<typename T> T* index_vector( size_t n ) {
		T *ivec = new T[ n ];
		T Tn = (T)n;
		for (T i = 0; i < Tn; i++) {
			ivec[ i ] = i;
		}
		return ivec;
	}

	// linearly spaced vector from firstval to lastval
	template<typename T>
	void linspace( T firstval, T lastval, size_t N,
			T *&vec ) {
		T stepsize = (lastval - firstval) / ((T)(N - 1));
		std::memset( vec, 0, N*sizeof(T) );
		for (size_t i = 0; i < N; i++) {
			vec[i] = firstval + ((T)i) * stepsize;
		}
	}
	template<typename T>
	T* linspace( T firstval, T lastval, size_t N ) {
		T *vec = NCPA::zeros<T>( N );
		NCPA::linspace<T>( firstval, lastval, N, vec );
		return vec;
	}

	// logarithmically spaced vector from a to b (NOT from 10^a to 10^b)
	template<typename T>
	void logspace( T a, T b, size_t k, T *&ls ) {
		T la = std::log10(a);
		T lb = std::log10(b);
		NCPA::linspace( la, lb, k, ls );
		for (size_t i = 0; i < k; i++) {
			ls[ i ] = std::pow( (T)(10.0), ls[ i ] );
		}
	}

	template<typename T>
	T* logspace( T firstval, T lastval, size_t N ) {
		T *vec = NCPA::zeros<T>( N );
		NCPA::logspace<T>( firstval, lastval, N, vec );
		return vec;
	}


	template<typename T> T** allocate_matrix( size_t nr, size_t nc ) {
		T **v;
		v = new T* [nr];
		for (size_t i = 0; i < nr; i++) {
			v[ i ] = new T[ nc ];
			std::memset( v[ i ], 0, nc * sizeof(T) );
		}
		return v;
	}

	template<typename T> void free_matrix( T **v, size_t nr, size_t nc ) {
		for (size_t i = 0; i < nr; i++) {
			delete [] v[ i ];
		}
		delete [] v;
	}

	template<typename T> void free_matrix_and_contents( T **v, size_t nr, size_t nc ) {
		for (size_t i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++) {
				delete v[ i ][ j ];
			}
			delete [] v[ i ];
		}
		delete [] v;
	}

	template<typename T> void fill_matrix( T **v, size_t nr, size_t nc, const T &val ) {
		for (size_t i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++) {
				v[ i ][ j ] = val;
			}
		}
	}

	template<typename T>
	size_t find_closest_index( T *z, size_t NZ, T zs ) {
		if (NZ == 1) {
			return 0;
		}
		double diff = 0.0, mindiff = DBL_MAX;
		size_t tmpind = 0;

		for (size_t i = 0; i < NZ; i++) {
			diff = std::abs( ((double)z[i]) - ((double)zs) );
			if (diff < mindiff) {
				tmpind = i;
				mindiff = diff;
			}
		}

		return tmpind;
	}

	template<typename T>
	size_t find_closest_index( std::vector<T> z, T zs ) {
		if (z.size() == 1) {
			return 0;
		}
		T diff = 0.0, mindiff = 9999999999999;
		size_t tmpind = 0, i;
		for (i = 0; i < z.size(); i++) {
			diff = std::abs( z[i] - zs );
			if (diff < mindiff) {
				tmpind = i;
				mindiff = diff;
			}
		}
		return tmpind;
	}

	template<typename T> T*** matrix3d( size_t nd1, size_t nd2, size_t nd3 ) {
		size_t i, j;
		T ***p = new T** [ nd1 ];

		for (i=0; i < nd1; i++) {
		  p[ i ] = new T*[ nd2 ];
		  for (j=0; j < nd2; j++) {
		  	p[ i ][ j ] = new T[ nd3 ];
		  	std::memset( p[ i ][ j ], 0, nd3 * sizeof(T) );
		  }
		}
		return p;
	}

	template<typename T> void free_matrix3d( T ***data, size_t nd1, size_t nd2, size_t nd3 ) {
		size_t i, j;
		for (i=0; i < nd1; ++i) {
		  if (data[i] != NULL) {
		      for (j=0; j < nd2; ++j) {
		      	  delete [] data[i][j];
		      }
		      delete [] data[i];
		  }
		}
		delete [] data;
	}


	template<typename T> void free_matrix3d_and_contents( T ***data,
				size_t nd1, size_t nd2, size_t nd3 ) {
		size_t i, j, k;
		for (i=0; i < nd1; ++i) {
		  if (data[i] != NULL) {
		      for (j=0; j < nd2; ++j) {
	      	  	  for (k = 0; i < nd3; k++) {
	      	  	  	  delete data[i][j][k];
	      	  	  }
		          delete [] data[i][j];
		      }
		      delete [] data[i];
		  }
		}
		delete [] data;
	}

	template<typename T> void fill_matrix3d( T ***v, size_t nr, size_t nc, size_t nz,
				const T &val ) {
		for (size_t i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++) {
				for (size_t k = 0; k < nz; k++) {
					v[ i ][ j ][ k ] = val;
				}
			}
		}
	}

	template<typename T, typename U> void print_2_columns( std::ostream &out, size_t nz,
		T *col1, U *col2, std::string del = " " ) {
		for (size_t i = 0; i < nz; i++) {
			out << col1[i] << del << col2[i] << std::endl;
		}
	}

	template<typename T>
	T trapz( size_t N, T *xvec, T *yvec ) {
		assert(N > 1);
		T sum = 0.0;
		for (size_t i = 1; i < N; i++) {
			sum += (yvec[i] + yvec[i-1]) * 0.5 * (xvec[i] - xvec[i-1]);
		}
		return sum;
	}

	template<typename T,typename U>
	void vector_scale( size_t N, U *in, T factor, U *&out ) {
		U *tempvec = NCPA::zeros<U>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = in[ i ] * factor;
		}
		std::memcpy( out, tempvec, N*sizeof(U) );
		delete [] tempvec;
	}

	template<typename T> void vector_multiply( size_t N, T *v1, T *v2,
		T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] * v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	template<typename T> void vector_add( size_t N, T *v1, T *v2,
		T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] + v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	template<typename T> void vector_subtract( size_t N, T *v1, T *v2,
		T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] - v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	template<typename T> void vector_divide( size_t N, T *v1, T *v2,
		T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] / v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	template<typename T> void circshift( T *in, size_t N, int shift,
			T *&out ) {
		assert(((size_t)abs(shift)) < N);
		size_t i;
		T *tempvec = NCPA::zeros<T>( N );
		if (shift < 0) {
			// move from the back to the front
			size_t negshift = (size_t)(-shift);
			for (i = negshift; i < N; i++) {
				tempvec[ i - negshift ] = in[ i ];
			}
			for (i = 0; i < negshift; i++) {
				tempvec[ N - negshift + i ] = in[ i ];
			}
		} else if (shift > 0) {
			// move from the front to the back
			for (i = 0; i < N-shift; i++) {
				tempvec[ i + shift ] = in[ i ];
			}
			for (i = N - shift; i < N; i++) {
				tempvec[ i - (N - shift) ] = in[ i ];
			}
		} else {
			std::memcpy( tempvec, in, N*sizeof(T) );
		}

		// copy over so it can be done in place if desired
		std::memcpy( out, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	template<typename T> void reverse( T *in, size_t N, T *&out ) {
		T *tempvec = NCPA::zeros<T>( N );
		size_t j;
		for (j = 0; j < N; j++) {
			tempvec[ N-1-j ] = in[ j ];
		}
		std::memcpy( out, tempvec, N * sizeof(T) );

		delete [] tempvec;
	}




	// Class for 3-D matrices with 1-D storage behind the scenes.  Cuts a lot
	// of looping out
	template<class T> class StorageOptimizedMatrix3D {
	protected:
		T *data;
		size_t dim1_, dim2_, dim3_;

		size_t indices_3_to_1( size_t ind1, size_t ind2, size_t ind3 ) const {
			// size_t ind1d = ind1 * dim2_ * dim3_ + ind2 * dim3_ + ind3;
			// std::cout << "[" << ind1 << "," << ind2 << "," << ind3 << "] = "
			// 		  << ind1d << std::endl;
			// return ind1d;
			return ind1 * dim2_ * dim3_ + ind2 * dim3_ + ind3;
		}

	public:
		StorageOptimizedMatrix3D( size_t dim1, size_t dim2, size_t dim3 ) {
			data = zeros<T>( dim1 * dim2 * dim3 );
			dim1_ = dim1;
			dim2_ = dim2;
			dim3_ = dim3;
		}

		~StorageOptimizedMatrix3D() {
			delete [] data;
		}

		void set( size_t ind1, size_t ind2, size_t ind3, T val ) {
			data[ indices_3_to_1( ind1, ind2, ind3 ) ] = val;
		}

		T get( size_t ind1, size_t ind2, size_t ind3 ) const {
			return data[ indices_3_to_1( ind1, ind2, ind3 ) ];
		}
	};
}






#endif
