/*
	UTIL.H: NCPA Utility functions library

	General:
	within(double,double,size_t): bool


	Array manipulation:
	circshift<T>(T*,size_t,int,t*&): void
	fill_matrix<T>(T**,size_t,size_t,const T&): void
	fill_matrix3d<T>(T***,size_t,size_t,size_t,const T&): void
	find_closest_index<T>(T*,size_t,T): size_t
	find_closest_index<T>(std::vector<T>,T): size_t


	Conversion:
	cart2pol(double,double,double&,double&): void
	void complex2double( size_t n, const std::complex<double> *in,
		double *real, double *imag );
	deg2rad<T>(T): T
	void double2complex( size_t n, const double *real, std::complex<double> *out );
	void double2complex( size_t n, const double *real, const double *imag,
		std::complex<double> *out );
	rad2deg<T>(T): T

	File interaction:
	count_rows_arbcol(const std::string&): int
	fexists(const char*): bool

	Math:
	evalpoly<T>(size_t,const T*,T): T
	evalpoly<T>(const std::vector<T>&,T): T
	linear_interp

	Memory allocation:
	allocate_matrix<T>(size_t,size_t): T**
	c3Darray(size_t,size_t,size_t): std::complex<double>*** DEPRECATED, use matrix3d<std::complex<double>>
	cmatrix(long,long): std::complex<double>** DEPRECATED, USE allocate_matrix<std::complex<double>>
	dmatrix(long,long): double** DEPRECATED, USE allocate_matrix<double>
	free_c3Darray(std::complex<double>***,size_t,size_t): void DEPRECATED, USE free_matrix3d<std::complex<double>>
	free_cmatrix(std::complex<double>**,long,long): int DEPRECATED, use free_matrix<std::complex<double>>
	free_dmatrix(double**,long,long): int DEPRECATED, use free_matrix<double>
	free_imatrix(int**,long,long): int DEPRECATED, use free_matrix<int>
	free_matrix<T>(T**,size_t,size_t): void
	free_matrix_and_contents<T>(T**,size_t,size_t): void
	free_matrix3d<T>(T***,size_t,size_t,size_t): void
	free_matrix3d_and_contents<T>(T***,size_t,size_t,size_t): void
	imatrix(long,long): int** DEPRECATED, USE allocate_matrix<int>
	index_vector<T>(size_t): T*
	index_vector<T>(size_t, T): T*



	String Manipulation:
	deblank(const std::string&): std::string
	deblank(const std::string&, const std::string&): std::string



*/




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
#include <cmath>

namespace NCPA {

	/**
	Dynamically allocates a new array and sets all elements to zero
	before returning it.
	@brief Returns a new array of all zeros.
	@param n The size of the array.
	@returns A pointer to the newly-allocated array.
	*/
	template<typename T>
	T *zeros( size_t n ) {
		T *out = new T[ n ]();
		return out;
	}

	/**
	@brief Dynamically allocates and returns a two-dimensional array.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns A pointer to the newly-allocated array.
	*/
	template<typename T>
	T** allocate_matrix( size_t nr, size_t nc ) {
		T **v;
		v = new T* [nr];
		for (size_t i = 0; i < nr; i++) {
			v[ i ] = new T[ nc ]();
			//std::memset( v[ i ], 0, nc * sizeof(T) );
		}
		return v;
	}

	/**
	Dynamically allocates a 3-D array of complex values and returns
	a pointer to the array.
	@brief Allocates and returns a std::complex<double>***.
	@param xlen The first dimension of the array.
	@param ylen The second dimension of the array.
	@param zlen The third dimension of the array.
	@returns A pointer to the first value of the array.
	\deprecated{Deprecated in favor of the NCPA::matrix3d<> template.}
	*/
	std::complex<double> ***c3Darray(size_t xlen, size_t ylen, size_t zlen);

	/**
	@brief Converts from Cartesian to polar coordinates.
	@param x The x Cartesian input coordinate.
	@param y The y Cartesian input coordinate.
	@param r The r polar output coordinate.
	@param theta_rad The theta output polar coordinate, in radians.
	*/
	void cart2pol( double x, double y, double &r, double &theta_rad );

	/**
	Circularly shifts the elements in an array X by K positions. If K is
    positive, then the values of X are circularly shifted from the
    beginning to the end. If K is negative, they are shifted from the
    end to the beginning.  The resultant shifted array is returned in
    a new dynamically-allocated array.
    @brief Circularly shifts array elements.
    @param X The array whose elements are to be shifted.
    @param N The size of the array.
    @param K The number of positions to shift.
    @param out The shifted array.  Can be the same as either input to perform in-place.
    */
	template<typename T>
	void circshift( T *X, size_t N, int K,
			T *out ) {
		assert(((size_t)abs(K)) < N);
		size_t i;
		T *tempvec = NCPA::zeros<T>( N );
		if (K < 0) {
			// move from the back to the front
			size_t negshift = (size_t)(-K);
			// first, move the last -K values to the front
			for (i = 0; i < negshift; i++) {
				tempvec[ i ] = X[ N-negshift+i ];
			}
			for (i = negshift; i < N; i++) {
				tempvec[ i ] = X[ i-negshift ];
			}
		} else if (K > 0) {
			// move from the front to the back
			// first, move the first K values to the back
			for (i = 0; i < K; i++) {
				tempvec[ N-K+i ] = X[ i ];
			}
			// Now, move the next N-K values to the front
			for (i = K; i < N; i++) {
				tempvec[ i-K ] = X[ i ];
			}
		} else {
			std::memcpy( tempvec, X, N*sizeof(T) );
		}

		// copy over so it can be done in place if desired
		std::memcpy( out, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}




	/**
	Dynamically allocates a 2-D array of complex values and returns
	a pointer to the array.
	@brief Allocates and returns a std::complex<double>**.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns A pointer to the first value of the array.
	\deprecated{Deprecated in favor of the NCPA::matrix<> class template
		(if matrix arithmetic capability is needed) or the
		NCPA::allocate_matrix<> template (if not).}
	*/
	std::complex<double> **cmatrix(long nr, long nc);

	/**
	 * Converts an array of complex<double>s into the equivalent arrays of doubles
	 * where one array represents the real parts and the other the imaginary parts.
	 * @brief Converts a complex<double> array to two double arrays
	 * @param n The size of the arrays
	 * @param in The complex array to decompose
	 * @param real The real parts
	 * @param imag The imaginary parts
	 */
	void complex2double( size_t n, const std::complex<double> *in,
			double *real, double *imag );

	/**
	@brief Counts the rows in a file.
	@input filename The filename to count the rows from.
	@returns The number of rows in the file.
	*/
	int count_rows_arbcol( const std::string& filename );

	/**
	Removes leading and trailing whitespace from a string.
	Whitespace is defined as spaces, tabs, or newline characters.
	@brief Removes leading and trailing whitespace from a string.
	@param orig The original string to deblank.
	@returns The original string with leading and trailing whitespace removed.
	*/
	std::string deblank( const std::string& orig );

	/**
	Removes leading and trailing whitespace from a string, using
	specified whitespace characters.
	@brief Removes leading and trailing whitespace from a string.
	@param orig The original string to deblank.
	@param whitespace A string containing all characters to treat as whitespace
	@returns The original string with leading and trailing whitespace removed.
	*/
	std::string deblank( const std::string& str,
		const std::string& whitespace );

	/**
	@brief Converts from degrees to radians.
	@param deg_in The input value in degrees.
	@returns The output value in radians.
	*/
	template<typename T>
	T deg2rad( T deg_in ) {
		return (T)(((double)deg_in) * M_PI / 180.0);
	}

	/**
	Dynamically allocates a 2-D array of double values and returns
	a pointer to the array.
	@brief Allocates and returns a double**.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns A pointer to the first value of the array.
	\deprecated{Deprecated in favor of the NCPA::matrix<> class template
		(if matrix arithmetic capability is needed) or the
		NCPA::allocate_matrix<> template (if not).}
	*/
	double **dmatrix(long nr, long nc);

	/**
	 * Converts an array of doubles into the equivalent array of complex<double>s
	 * with the imaginary value set to zero.
	 * @brief Converts a double array to a complex<double> array
	 * @param n The size of the arrays
	 * @param real The double array to convert
	 * @param out The complex<double> array to return
	 */
	void double2complex( size_t n, const double *real, std::complex<double> *out );

	/**
	 * Converts two arrays of doubles into the equivalent array of complex<double>s
	 * where one array represents the real parts and the other the imaginary parts.
	 * @brief Converts a double array to a complex<double> array
	 * @param n The size of the arrays
	 * @param real The double array to convert as the real parts
	 * @param imag The double array to convert as the real parts
	 * @param out The complex<double> array to return
	 */
	void double2complex( size_t n, const double *real, const double *imag,
			std::complex<double> *out );

	/**
	@brief Evaluate a polynomial from its coefficients.
	@param N the number of coefficients (not its degree).
	@param coeffs The coefficients in increasing order starting with power 0 (i.e. the constant).
	@param x The value at which to evaluate the polynomial.
	@returns The calculated value of the polynomial at x.
	*/
	template<typename T>
	T evalpoly( size_t N, const T *coeffs, T x ) {
		T val = 0.0;
		for (size_t i = 0; i < N; i++) {
			val += coeffs[i] * std::pow( x, (double)i );
		}
		return val;
	}

	/**
	@brief Evaluate a polynomial from its coefficients.
	@param coeffs The coefficients in increasing order starting with power 0 (i.e. the constant).
	@param x The value at which to evaluate the polynomial.
	@returns The calculated value of the polynomial at x.
	*/
	template<typename T>
	T evalpoly( const std::vector<T> &coeffs, T x ) {
		T val = 0.0;
		for (size_t i = 0; i < coeffs.size(); i++) {
			val += coeffs[i] * std::pow( x, (double)i );
		}
		return val;
	}

	/**
	Determines if a given filename represents a readable file by
	attempting to open it and determining if it is in a good state.
	@brief Determines if a given filename represents a readable file.
	@param filename The name of the file to attempt to open.
	@returns true if the file can be read on open, false otherwise.
	*/
	bool fexists( const char *filename );

	/**
	@brief Sets each element of a 2-D array to a constant value.
	@param v The array to fill.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@param val The value to set each element to.
	*/
	template<typename T>
	void fill_matrix( T **v, size_t nr, size_t nc, const T &val ) {
		for (size_t i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++) {
				v[ i ][ j ] = val;
			}
		}
	}

	/**
	@brief Sets each element of a three-dimensional array to a constant value.
	@param v The array to fill.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@param nz The third dimension of the array.
	@param val The value to set each element to.
	*/
	template<typename T>
	void fill_matrix3d( T ***v, size_t nr, size_t nc, size_t nz,
				const T &val ) {
		for (size_t i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++) {
				for (size_t k = 0; k < nz; k++) {
					v[ i ][ j ][ k ] = val;
				}
			}
		}
	}

	/**
	 * Finds and returns the indices of the array elements before and after the supplied
	 * value.
	 */
	template<typename T>
	bool find_interval_inclusive( T* z, size_t NZ, T val, size_t &bottom, size_t &top ) {
		double *it = std::lower_bound(z, z+NZ, val);
		size_t diff = it - z;
		if (diff == 0) {
			bottom = 0;
			if (val < z[0]) {
				top = 0;
				return false;
			} else {
				top = 1;
				return true;
			}
		} else if (diff == NZ) {
			top = NZ;
			if (val > z[NZ-1]) {
				bottom = NZ;
				return false;
			} else {
				bottom = NZ-1;
				return true;
			}
		} else {
			bottom = diff - 1;
			top = diff;
			return true;
		}
	}

	/**
	Finds and returns the index of the array element that is closest
	in value to a supplied reference value.  Minimizes abs( z - zref ).
	@brief Finds the index of the closest array element to a value.
	@param z The array to check.
	@param NZ The size of the array.
	@param zs The reference value.
	@returns The index of the array element closest to the reference value.
	*/
	template<typename T>
	size_t find_closest_index( T *z, size_t NZ, T zs ) {
		if (NZ == 1) {
			return 0;
		}
		double diff = 0.0, mindiff = DBL_MAX;
		size_t tmpind = 0;

		for (size_t i = 0; i < NZ; i++) {
			diff = std::fabs( ((double)z[i]) - ((double)zs) );
			if (diff < mindiff) {
				tmpind = i;
				mindiff = diff;
			}
		}

		return tmpind;
	}

	/**
	Finds and returns the index of the vector element that is closest
	in value to a supplied reference value.  Minimizes abs( z - zref ).
	@brief Finds the index of the closest vector element to a value.
	@param z The vector to check.
	@param zs The reference value.
	@returns The index of the vector element closest to the reference value.
	*/
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


	/**
	Frees a dynamically-allocated 3-D array of complex values.
	@brief Frees a std::complex<double>***.
	@param data The pointer to free.
	@param xlen The first dimension of the array.
	@param ylen The second dimension of the array.
	\deprecated{Deprecated in favor of NCPA::free_matrix3d<> template.}
	*/
	void free_c3Darray(std::complex<double> ***data, size_t xlen,
		size_t ylen);

	/**
	Frees a dynamically-allocated 2-D array of complex values.
	@brief Frees a std::complex<double>**.
	@param v The pointer to free.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns 0 on success.
	\deprecated{Deprecated in favor of NCPA::free_matrix<> template.}
	*/
	int free_cmatrix(std::complex<double> **v, long nr, long nc);

	/**
	Frees a dynamically-allocated 2-D array of double values.
	@brief Frees a double**.
	@param v The pointer to free.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns 0 on success.
	\deprecated{Deprecated in favor of NCPA::free_matrix<> template.}
	*/
	int free_dmatrix(double** v, long nr, long nc);

	/**
	Frees a dynamically-allocated 2-D array of int values.
	@brief Frees an int**.
	@param v The pointer to free.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns 0 on success.
	\deprecated{Deprecated in favor of NCPA::free_matrix<> template.}
	*/
	int free_imatrix(int** v, long nr, long nc);

	/**
	@brief Frees a dynamically allocated 2-D array.
	@param v The array to free.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	*/
	template<typename T>
	void free_matrix( T **v, size_t nr, size_t nc ) {
		for (size_t i = 0; i < nr; i++) {
			delete [] v[ i ];
		}
		delete [] v;
	}

	/**
	Frees a dynamically-allocated 2-D array of objects, and calls
	delete on each element in the process.
	@brief Frees a dynamically allocated 2-D array and its contents.
	@param v The array to free.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	*/
	template<typename T>
	void free_matrix_and_contents( T **v, size_t nr, size_t nc ) {
		for (size_t i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++) {
				delete v[ i ][ j ];
			}
			delete [] v[ i ];
		}
		delete [] v;
	}

	/**
	@brief Frees a dynamically-allocated three-dimensional array.
	@param data The array to free.
	@param nd1 The first dimension of the array.
	@param nd2 The second dimension of the array.
	@param nd3 The third dimension of the array.
	*/
	template<typename T>
	void free_matrix3d( T ***data, size_t nd1, size_t nd2, size_t nd3 ) {
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

	/**
	Frees a dynamically-allocated three-dimensional array, and calls
	delete on each element as it does so.
	@brief Frees a dynamically-allocated three-dimensional array and its contents.
	@param data The array to free.
	@param nd1 The first dimension of the array.
	@param nd2 The second dimension of the array.
	@param nd3 The third dimension of the array.
	*/
	template<typename T>
	void free_matrix3d_and_contents( T ***data,
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

	/**
	Dynamically allocates a 2-D array of int values and returns
	a pointer to the array.
	@brief Allocates and returns an int**.
	@param nr The first dimension of the array.
	@param nc The second dimension of the array.
	@returns A pointer to the first value of the array.
	\deprecated{Deprecated in favor of the NCPA::matrix<> class template
		(if matrix arithmetic capability is needed) or the
		NCPA::allocate_matrix<> template (if not).}
	*/
	int **imatrix(long nr, long nc);

	/**
	Dynamically allocates a new array and sets the elements to
	(0, 1, 2, ..., n-1).
	@brief Returns a new array of index values.
	@param n The size of the array.
	@returns A pointer to the newly-allocated array.
	*/
	template<typename T>
	T *index_vector( size_t n ) {
		T *ivec = new T[ n ];
		for (size_t i = 0; i < n; i++) {
			ivec[ i ] = (T)i;
		}
		return ivec;
	}

	/**
	Dynamically allocates a new array and sets the elements to
	(a, a+1, a+2, ..., a+n-1).
	@brief Returns a new array of index values.
	@param n The size of the array.
	@param a The offset from zero to use for each value.
	@returns A pointer to the newly-allocated array.
	*/
	template<typename T>
	T *index_vector( size_t n, T a ) {
		T *ivec = new T[ n ];
		for (size_t i = 0; i < n; i++) {
			ivec[ i ] = (T)i + a;
		}
		return ivec;
	}


	/**
	Performs a simple linear interpolation.  Returns the value at x
	between (x1,y1) and (x2,y2).
	@brief Performs a simple linear interpolation.
	@param x1 The x coordinate of the start point.
	@param y1 The y coordinate of the start point.
	@param x2 The x coordinate of the end point.
	@param y2 The y coordinate of the end point.
	@param x The x value at which to perform the interpolation.
	@returns The value of y at x.
	*/
	double linearInterp( double x1, double y1, double x2, double y2, double x );

	/**
	Returns an array of N values, linearly spaced between two
	supplied end points, placed in a pre-allocated array.
	@brief Populates a linearly spaced array of values.
	@param firstval The starting value of the array.
	@param lastval The ending value of the array.
	@param N The number of values to return.
	@param vec The array to place the values in.
	*/
	template<typename T>
	void linspace( T firstval, T lastval, size_t N,
			T *&vec ) {
		T stepsize = (lastval - firstval) / ((T)(N - 1));
		std::memset( vec, 0, N*sizeof(T) );
		for (size_t i = 0; i < N; i++) {
			vec[i] = firstval + ((T)i) * stepsize;
		}
	}

	/**
	Dynamically allocates an array of N values, linearly spaced
	between two supplied end points, and returns the new array.
	@brief Returns a new linearly spaced array of values.
	@param firstval The starting value of the array.
	@param lastval The ending value of the array.
	@param N The number of values to return.
	@returns The new array containing the values.
	*/
	template<typename T>
	T* linspace( T firstval, T lastval, size_t N ) {
		T *vec = NCPA::zeros<T>( N );
		NCPA::linspace<T>( firstval, lastval, N, vec );
		return vec;
	}

	// logarithmically spaced vector from a to b (NOT from 10^a to 10^b)
	/**
	Returns an array of N values, logarithmically spaced between two
	supplied end points a and b, placed in a pre-allocated array.  Note
	that this returns the values, not the exponents of the values as in
	some other logspace functions; in other words, it returns the array
	(a, ..., b) and NOT (10^a, ..., 10^b).
	@brief Populates a linearly spaced array of values.
	@param a The starting value of the array.
	@param b The ending value of the array.
	@param k The number of values to return.
	@param ls The array to place the values in.
	*/
	template<typename T>
	void logspace( T a, T b, size_t k, T *&ls ) {
		T la = std::log10(a);
		T lb = std::log10(b);
		NCPA::linspace( la, lb, k, ls );
		for (size_t i = 0; i < k; i++) {
			ls[ i ] = std::pow( (T)(10.0), ls[ i ] );
		}
	}

	/**
	Dynamically allocates and returns an array of N values, logarithmically
	spaced between two supplied end points a and b, placed in a
	pre-allocated array.
	Note that this returns the values, not the exponents of the values
	as in some other logspace functions; in other words, it returns
	the array (a, ..., b) and NOT (10^a, ..., 10^b).
	@brief Populates a linearly spaced array of values.
	@param a The starting value of the array.
	@param b The ending value of the array.
	@param N The number of values to return.
	@returns The new array containing the values.
	*/
	template<typename T>
	T* logspace( T a, T b, size_t N ) {
		T *vec = NCPA::zeros<T>( N );
		NCPA::logspace<T>( a, b, N, vec );
		return vec;
	}

	/**
	@brief Dynamically allocates and returns a three-dimensional array.
	@param nd1 The first dimension of the array.
	@param nd2 The second dimension of the array.
	@param nd3 The third dimension of the array.
	@returns A pointer to the newly-allocated array.
	*/
	template<typename T>
	T*** matrix3d( size_t nd1, size_t nd2, size_t nd3 ) {
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

	/**
	@brief Returns the greater of two values.
	@param a First value.
	@param b Second value.
	@returns a if a > b, b otherwise.
	*/
	template<typename T>
	T max( T a, T b ) {
		return a > b ? a : b;
	}

	/**
	@brief Returns the minimum value from a vector.
	@param vals The vector holding the values to check.
	@returns The minimum value found in the vector.
	*/
	template<typename T>
	T max( std::vector< T > vals ) {
		T maxval = vals.front();
		for (typename std::vector<T>::const_iterator cit = vals.cbegin();
				cit != vals.cend(); ++cit) {
			maxval = NCPA::max( *cit, maxval );
		}
		return maxval;
	}

	/**
	@brief Returns the maximum value from an array.
	@param vals The array holding the values to check.
	@param size The size of the array.
	@returns The maximum value found in the array.
	*/
	template<typename T>
	T max( const T *vals, size_t size ) {
		T maxval = vals[ 0 ];
		for (size_t i = 1; i < size; i++) {
			maxval = NCPA::max( vals[i], maxval );
		}
		return maxval;
	}

	/**
	Computes the mean of a double array.
	@brief Computes the mean of a double array.
	@param d The array of doubles to average.
	@param n The size of the array.
	@returns The computed mean.
	*/
	double mean( double* d, size_t n );

	/**
	@brief Returns the lesser of two values.
	@param a First value.
	@param b Second value.
	@returns a if a < b, b otherwise.
	*/
	template<typename T>
	T min( T a, T b ) {
		return a < b ? a : b;
	}

	/**
	@brief Returns the minimum value from an array.
	@param vals The array holding the values to check.
	@param size The size of the array.
	@returns The minimum value found in the array.
	*/
	template<typename T>
	T min( const T *vals, size_t size ) {
		T minval = vals[ 0 ];
		for (size_t i = 1; i < size; i++) {
			minval = NCPA::min( vals[ i ], minval );
		}
		return minval;
	}



	/**
	@brief Returns the maximum value from a vector.
	@param vals The vector holding the values to check.
	@returns The maximum value found in the vector.
	*/
	template<typename T>
	T min( std::vector< T > vals ) {
		T minval = vals.front();
		for (typename std::vector<T>::const_iterator cit = vals.cbegin();
				cit != vals.cend(); ++cit) {
			minval = NCPA::min( *cit, minval );
		}
		return minval;
	}

	/**
	Shifts a point in Cartesian coordinates to a new position
	relative to a different origin.
	@brief Returns coordinates relative to a new origin.
	@param old_x The x coordinate relative to (0,0).
	@param old_y The y coordinate relative to (0,0).
	@param x_new_origin The x coordinate of the new origin.
	@param y_new_origin The y coordinate of the new origin.
	@param new_x The new x coordinate.
	@param new_y The new y coordiante.
	*/
	template<typename T>
	void move_origin( T old_x, T old_y,
			T x_new_origin, T y_new_origin, T &new_x, T &new_y ) {
		new_x = old_x - x_new_origin;
		new_y = old_y - y_new_origin;
	}

	/**
	Shifts an array of Cartesian points from an origin at (0,0) to
	a new origin, and returns the coordinates relative to the new
	origin.
	@brief Returns coordinates relative to a new origin.
	@param npts The number of points to move.
	@param old_x The x coordinates relative to (0,0)
	@param old_y The y coordinates relative to (0,0)
	@param x_new_origin The x coordinate of the new origin
	@param y_new_origin The y coordinate of the new origin
	@param new_x The array to hold the new x coordinates.
	@param new_y The array to hold the new y coordinates.
	*/
	template<typename T>
	void move_origin( size_t npts, const T *old_x, const T *old_y,
			T x_new_origin, T y_new_origin, T *new_x, T *new_y ) {
		for (size_t i = 0; i < npts; i++) {
			move_origin<T>( old_x[ i ], old_y[ i ],
				x_new_origin, y_new_origin, new_x[ i ], new_y[ i ] );
		}
	}

	/**
	For input v, returns the lowest positive integer n such that
	2^n >= v.
	@brief Returns the next integer power of 2.
	@param v The input number to start at.
	@returns The next integer power of 2.
	*/
	size_t nextpow2( size_t v );

	/**
	@brief Normalizes an azimuth into the range [0,360)
	@param in The azimuth value to normalize, in degrees.
	@returns The normalized azimuth, in degrees.
	*/
	template<typename T>
	T normalizeAzimuth( T in ) {
		double out = (double)in;
		while (out < 0.0) {
			out += 360.0;
		}
		while (out >= 360.0) {
			out -= 360.0;
		}
		return (T)out;
	}

	/**
	@brief Normalizes an azimuth into the range [0,2*pi)
	@param in The azimuth value to normalize, in radians.
	@returns The normalized azimuth, in radians.
	*/
	template<typename T>
	T normalizeAzimuthRadians( T in ) {
		double out = (double)in;
		while (out < 0.0) {
			out += 2.0 * M_PI;
		}
		while (out >= 2.0*M_PI) {
			out -= 2.0 * M_PI;
		}
		return (T)out;
	}

	/**
	@brief Converts from polar to Cartesian coordinates.
	@param r The r polar input coordinate.
	@param theta_rad The theta input polar coordinate, in radians.
	@param x The x Cartesian output coordinate.
	@param y The y Cartesian output coordinate.
	*/
	void pol2cart( double r, double theta_rad, double &x, double &y );

	/**
	@brief Prints two arbitrary columns to a stream with specified delimiter.
	@param out The stream to output to.
	@param nz The number of values in each array.
	@param col1 The first column.
	@param col2 The second column.
	@param del The delimiter.  Defaults to a single space.
	*/
	template<typename T, typename U>
	void print_2_columns( std::ostream &out, size_t nz,
		T *col1, U *col2, std::string del = " " ) {
		for (size_t i = 0; i < nz; i++) {
			out << col1[i] << del << col2[i] << std::endl;
		}
	}

	/**
	@brief Converts from radians to degrees.
	@param rad_in The input value in radians.
	@returns The output value in degrees.
	*/
	template<typename T>
	T rad2deg( T rad_in ) {
		return (T)(((double)rad_in) * 180.0 / M_PI);
	}

	/**
	Generates a vector of N random numbers in the range [0,scale).
	@brief Generates a vector of random numbers.
	@param N Number of random values to generate.
	@param scale The upper end of the range.  Defaults to 1.0.
	@returns A vector of random numbers.
	*/
	std::vector<double> random_numbers( size_t N, double scale = 1.0 );

	/**
	@brief Reads a matrix from a text file, one row per line.
	@param filename The file to read from.
	@param contents A dynamically-allocated double[ncols][nrows].
	@param ncols The number of columns read from the file.
	@param nrows The number of rows read from the file.
	@param delimiters Delimiters to use between values.  Defaults to whitespace.
	*/
	void read_matrix_from_file( const std::string &filename,
		double **&contents, size_t &ncols, size_t &nrows,
		std::string delimiters = " \t\r\n" );

	/**
	@brief Reads columns from a text file.
	@param filename The file to read from.
	@param contents The columns from the file, one vector of strings per column.
	@param delimiters Delimiters to use between values.  Defaults to whitespace.
	*/
	void read_text_columns_from_file( const std::string &filename,
		std::vector< std::vector< std::string > > &contents,
		std::string delimiters = " \t\r\n" );


	void read_text_columns_from_file_with_header(
		const std::string &filename,
		std::vector< std::vector< std::string > > &contents,
		std::vector< std::string > &headerlines,
		const std::string &delimiters = " \t\r\n",
		const std::string &headerchars = "#" );

	/**
	@brief Reverses the order of an array.
	@param in The input array.
	@param N The size of the array.
	@param out Holds the output array. Can be the same as the input.
	*/
	template<typename T>
	void reverse( T *in, size_t N, T *&out ) {
		T *tempvec = NCPA::zeros<T>( N );
		size_t j;
		for (j = 0; j < N; j++) {
			tempvec[ N-1-j ] = in[ j ];
		}
		std::memcpy( out, tempvec, N * sizeof(T) );

		delete [] tempvec;
	}


	/**
	Acts as std::getline, but checks for all variations of newline
	characters \r, \n, and \r\n.
	@brief Gets a line, checking for several variants of newline.
	@param is The stream to read the line from
	@param s The line returned.
	@returns The original stream after reading the line.
	*/
	std::istream &safe_getline( std::istream &is, std::string &s );

	/**
	Dynamically allocates a new array and sets all elements to the same
	value before returning it.
	@brief Returns a new array of a single value.
	@param n The size of the array.
	@param val The value to set all elements to.
	@returns A pointer to the newly-allocated array.
	*/
	template<typename T>
	T *single_valued_vector( size_t n, T val ) {
		T *out = new T[ n ];
		for (size_t i = 0; i < n; i++) {
			out[ i ] = val;
		}
		return out;
	}

	/**
	Splits a string into tokens using supplied delimiter characters.
	Similar to the Perl split() function.
	@brief Breaks a string into tokens using delimiter characters.
	@param input The string to split.
	@param delimiters The delimiter characters.  Defaults to whitespace.
	@returns A vector of the string tokens.
	*/
	std::vector< std::string > split( std::string input,
		std::string delimiters = " \t\n\r" );

	// Class for 3-D matrices with 1-D storage behind the scenes.  Cuts a lot
	// of looping out
	/**
	A three-dimensional array wrapper with optimized storage.
	A class that represents a three-dimensional array, but maintains
	internal storage as a one-dimensional array with array index
	translation handled internally, so that initialization and cleanup
	are not looped.
	*/
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

		/**
		@brief Basic constructor.
		@param dim1 First logical dimension of array.
		@param dim2 Second logical dimension of array.
		@param dim3 Third logical dimension of array.
		*/
		StorageOptimizedMatrix3D( size_t dim1, size_t dim2, size_t dim3 ) {
			data = zeros<T>( dim1 * dim2 * dim3 );
			dim1_ = dim1;
			dim2_ = dim2;
			dim3_ = dim3;
		}

		/**
		Basic destructor.
		*/
		~StorageOptimizedMatrix3D() {
			delete [] data;
		}

		/**
		@brief Set a value of the array.
		@param ind1 First index.
		@param ind2 Second index.
		@param ind3 Third index.
		@param val The value to set.
		*/
		void set( size_t ind1, size_t ind2, size_t ind3, T val ) {
			data[ indices_3_to_1( ind1, ind2, ind3 ) ] = val;
		}

		/**
		@brief Return a value of the array.
		@param ind1 First index.
		@param ind2 Second index.
		@param ind3 Third index.
		@returns The value of the array at (ind1,ind2,ind3).
		*/
		T get( size_t ind1, size_t ind2, size_t ind3 ) const {
			return data[ indices_3_to_1( ind1, ind2, ind3 ) ];
		}
	};

	/**
	Converts a numerical time to its equivalent time string in the
	format yyyy-mm-dd hh:mm:ss.mmm
	@brief Returns a human-readable time from an epoch time.
	@param d The time to print, as an epoch time_t value.
	@returns The string representation of the time.
	*/
	std::string timeAsString( double d );

	/**
	Converts the supplied string to lower case.
	@brief Converts to lower case.
	@param in The string to convert.
	@returns The input string converted into all lower case.
	*/
	void toLowerCase( char *in );

	/**
	Converts the supplied string to lower case.
	@brief Converts to lower case.
	@param in The string to convert.
	@returns The input string converted into all lower case.
	*/
	std::string toLowerCase( const std::string in );

	/**
	Converts the supplied string to upper case.
	@brief Converts to upper case.
	@param in The string to convert.
	@returns The input string converted into all upper case.
	*/
	void toUpperCase( char *in );

	/**
	Converts the supplied string to upper case.
	@brief Converts to upper case.
	@param in The string to convert.
	@returns The input string converted into all upper case.
	*/
	std::string toUpperCase( const std::string in );

	/**
	@brief Peforms a simple trapezoidal integration.
	@param N The number of points to integrate over.
	@param xvec The x values.
	@param yvec The y values.
	@returns The computed integral.
	*/
	template<typename T>
	T trapz( size_t N, T *xvec, T *yvec ) {
		assert(N > 1);
		T sum = 0.0;
		for (size_t i = 1; i < N; i++) {
			sum += (yvec[i] + yvec[i-1]) * 0.5 * (xvec[i] - xvec[i-1]);
		}
		return sum;
	}

	/**
	Adds two arrays together element-wise, returning the
	sum in a supplied array.
	@brief Performs element-wise array addition.
	@param N The number of points in the array.
	@param v1 The first array to add.
	@param v2 The second array to add.
	@param v12 The array to hold the sum.  Can be the same as either input array, in which case the values are replaced.
	*/
	template<typename T>
	void vector_add( size_t N, T *v1, T *v2, T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] + v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	/**
	Divides one array by another element-wise, returning the
	quotient in a supplied array.
	@brief Performs element-wise array division.
	@param N The number of points in the array.
	@param v1 The array to divide.
	@param v2 The array to divide by.
	@param v12 The array to hold the quotient.  Can be the same as either input array, in which case the values are replaced.
	*/
	template<typename T>
	void vector_divide( size_t N, T *v1, T *v2,
		T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] / v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	/**
	Multiplies two arrays together element-wise, returning the
	product in a supplied array.
	@brief Performs element-wise array multiplication.
	@param N The number of points in the array.
	@param v1 The first array to multiply.
	@param v2 The second array to multiply.
	@param v12 The array to hold the product.  Can be the same as either input array, in which case the values are replaced.
	*/
	template<typename T>
	void vector_multiply( size_t N, T *v1, T *v2, T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] * v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}


	/**
	Scales an array of values by a constant value, returning the
	result in a dynamically-allocated array.
	@brief Scales an array by a constant value.
	@param N The number of points in the array.
	@param in The array to scale.
	@param factor The factor to scale by.
	@param out The new, dynamically-allocated scaled array.
	*/
	template<typename T,typename U>
	void vector_scale( size_t N, U *in, T factor, U *&out ) {
		U *tempvec = NCPA::zeros<U>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = in[ i ] * factor;
		}
		std::memcpy( out, tempvec, N*sizeof(U) );
		delete [] tempvec;
	}

	/**
	Scales an array of values in-place by a constant value
	@brief Scales an array by a constant value in place.
	@param N The number of points in the array.
	@param in The array to scale.
	@param factor The factor to scale by.
	*/
	template<typename T,typename U>
	void vector_scale( size_t N, U *in, T factor ) {
		for (size_t i = 0; i < N; i++) {
			in[ i ] *= factor;
		}
	}

	/**
	Subtracts one array from another element-wise, returning the
	difference in a supplied array.
	@brief Performs element-wise array subtraction.
	@param N The number of points in the array.
	@param v1 The array to subtract from.
	@param v2 The array to subtract.
	@param v12 The array to hold the difference.  Can be the same as either input array, in which case the values are replaced.
	*/
	template<typename T>
	void vector_subtract( size_t N, T *v1, T *v2, T *&v12 ) {
		T *tempvec = NCPA::zeros<T>(N);
		for (size_t i = 0; i < N; i++) {
			tempvec[ i ] = v1[i] - v2[i];
		}
		std::memcpy( v12, tempvec, N*sizeof(T) );
		delete [] tempvec;
	}

	/**
	@brief Compares two doubles to a given number of decimal places
	@param val1 First value
	@param val2 Second value
	@param precision Decimal places to compare
	@returns true if the numbers are equal to that many decimal places, false otherwise
	 */
	bool within( double val1, double val2, size_t precision );



}






#endif
