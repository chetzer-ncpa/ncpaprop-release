#include <cstdio>
#include <cmath>
#include <ctime>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <complex>
#include <cctype>
#include <cstring>
#include <vector>
#include <sstream>
#include <cfloat>
#include <random>
#include "util.h"

#ifndef PI
#define PI 3.141592653589793
#endif


std::string NCPA::timeAsString(double d) {
    time_t temptime = (time_t)d;
    tm* uttime = std::gmtime( &temptime );
    double ipart, fpart;
    fpart = std::modf( d, &ipart);

    char* holder = new char[ 28 ];
    std::sprintf(holder,"%4d-%02d-%02d %02d:%02d:%02d.%03d GMT",
            uttime->tm_year+1900, uttime->tm_mon+1, uttime->tm_mday,
            uttime->tm_hour, uttime->tm_min, uttime->tm_sec,
            (int)(round(fpart * 1000)) );
    std::string s = holder;
    return s;
}

bool NCPA::fexists( const char *filename ) {
	std::ifstream ifile( filename );
	bool tf = ifile.good();
	ifile.close();
	return tf;
}

// Acts like std::getline, but checks for all three permutations of EOL characters
std::istream &NCPA::safe_getline( std::istream &is, std::string &s ) {

	char ch;
	s.clear();

	// Keep going until you get to a linefeed character
	while (is.get(ch) && ch != '\n' && ch != '\r')
		s += ch;

	// DOS systems use 2 consecutive characters, so make sure there isn't another one lurking
	if (ch == '\r') {
		ch = is.peek();
		if (ch == '\n') {
			is.get(ch);
		}
	}
	return is;
}

// perform a simple linear interpolation
double NCPA::linearInterp( double x1, double y1, double x2, double y2, double x ) {
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

void NCPA::toLowerCase( char *in ) {
	unsigned int len = std::strlen(in);
	for (unsigned int i = 0; i < len; i++) {
		in[ i ] = (char)std::tolower( in[ i ] );
	}
}

void NCPA::toUpperCase( char *in ) {
	unsigned int len = std::strlen(in);
	for (unsigned int i = 0; i < len; i++) {
		in[ i ] = (char)std::toupper( in[ i ] );
	}
}

std::string NCPA::toLowerCase( const std::string in ) {
	std::ostringstream oss("");
	for (unsigned int i = 0; i < in.length(); i++) {
		oss << (char)std::tolower( in[ i ] );
	}
	return oss.str();
}

std::string NCPA::toUpperCase( const std::string in ) {
	std::ostringstream oss("");
	for (unsigned int i = 0; i < in.length(); i++) {
		oss << (char)std::toupper( in[ i ] );
	}
	return oss.str();
}

std::vector< std::string > NCPA::split( std::string input, std::string delimiters ) {
	// // lean on strtok for this functionality
	// char *holder = new char[ input.size() + 1 ];
	// std::memset( holder, '\0', input.size() + 1 );
	// std::strcpy( holder, input.c_str() );
	// char *tmp = strtok( holder, delimiters.c_str() );
	// std::vector< std::string > tokens;
	// tokens.clear();
	// while (tmp != NULL) {
	// 	std::string s( tmp );
	// 	tokens.push_back( s );
	// 	tmp = strtok(NULL,delimiters.c_str());
	// }
	// delete [] holder;
	// return tokens;

	std::vector< std::string > tokens;
	tokens.clear();
	size_t firstind, lastind;
	firstind = input.find_first_not_of( delimiters );
	lastind = input.find_first_of( delimiters, firstind );
	while( firstind != std::string::npos ) {
		tokens.push_back( input.substr( firstind, lastind - firstind ) );
		firstind = input.find_first_not_of( delimiters, lastind );
		lastind = input.find_first_of( delimiters, firstind );
	}
	return tokens;

}

// bool NCPA::checkAzimuthLimits( double toCheck, double target, double tolerance ) {

// 	return ( std::cos(target - toCheck) >= std::cos( tolerance ) );
// }

// double NCPA::normalizeAzimuth( double in ) {
// 	double out = in;
// 	while (out < 0) {
// 		out += 360;
// 	}
// 	while (out >= 360) {
// 		out -= 360;
// 	}
// 	return out;
// }

int** NCPA::imatrix( long nr, long nc ) {
	// allocate an int matrix
	int **v;
	v = new int* [nr];
	for (long i = 0; i < nr; i++) {
		v[ i ] = new int[ nc ];
	}
	return v;
}

int NCPA::free_imatrix( int **v, long nr, long nc ) {
	for (long i = 0; i < nr; i++) {
		delete [] v[ i ];
	}
	delete [] v;
	return 0;
}

// Following functions came from various <module>_lib files
double** NCPA::dmatrix(long nr, long nc) {
  // allocate a double matrix
  double **v;
  v = new double* [nr];
  for (long i=0; i<nr; i++) {
      v[i] = new double [nc];
  }
  return v;
}


int NCPA::free_dmatrix(double **v, long nr, long nc) {
  // free a double matrix
  for (long i=0; i<nr; i++) {
      delete [] v[i];
  }
  delete [] v;
  return 0;
}

std::complex<double> **NCPA::cmatrix(long nr, long nc) {
  // allocate a complex<double> matrix
  std::complex<double> **v;
  v = new std::complex<double>* [nr];
  for (long i=0; i<nr; i++) {
      v[i] = new std::complex<double> [nc];
  }
  return v;
}

int NCPA::free_cmatrix( std::complex<double> **v, long nr, long nc) {
  // free a complex<double> matrix
  for (long i=0; i<nr; i++) {
      delete [] v[i];
  }
  delete [] v;
  return 0;
}

std::complex<double> ***NCPA::c3Darray(size_t xlen, size_t ylen, size_t zlen)
{
  std::complex<double> ***p;
  size_t i, j;

  if ((p = (std::complex<double> ***) calloc(xlen, (sizeof(*p)))) == NULL) {
      perror("calloc 1");
      return NULL;
  }

  for (i=0; i < xlen; ++i)
      p[i] = NULL;

  for (i=0; i < xlen; ++i)
      if ((p[i] = (std::complex<double> **) calloc(ylen, (sizeof (*p[i])))) == NULL) {
          perror("calloc 2");
          NCPA::free_c3Darray(p, xlen, ylen);
          return NULL;
      }

  for (i=0; i < xlen; ++i)
      for (j=0; j < ylen; ++j)
          p[i][j] = NULL;

  for (i=0; i < xlen; ++i)
      for (j=0; j < ylen; ++j)
          if ((p[i][j] = (std::complex<double> *) calloc(zlen, (sizeof (*p[i][j])))) == NULL) {
              perror("calloc 3");
              NCPA::free_c3Darray(p, xlen, ylen);
              return NULL;
          }
  return p;
}


void NCPA::free_c3Darray(std::complex<double> ***data, size_t xlen, size_t ylen)
{
  size_t i, j;
  for (i=0; i < xlen; ++i) {
      if (data[i] != NULL) {
          for (j=0; j < ylen; ++j)
              free(data[i][j]);
          free(data[i]);
      }
  }
  free(data);
}

std::string NCPA::deblank( const std::string& str ) {
	return NCPA::deblank( str, " \t\n\r" );
}

std::string NCPA::deblank( const std::string& str, const std::string& whitespace ) {
	const size_t strBegin = str.find_first_not_of( whitespace );
	if (strBegin == std::string::npos) {
		return "";
	}

	const size_t strEnd = str.find_last_not_of( whitespace );
	const size_t strRange = strEnd - strBegin + 1;

	return str.substr( strBegin, strRange );
}

int NCPA::count_rows_arbcol(const std::string& filename) {
  int answer,c;
  FILE *f=fopen(filename.c_str(),"r");

  if(f==NULL) {
      std::ostringstream es;
      es << "file << " << filename << " could not be opened.\n";
      throw std::invalid_argument(es.str());
  }
  answer=0;
  //read_header(f);
  while((c=getc(f))!=EOF){
      if(c=='\n') answer=answer+1;
  }
  fclose(f);
  return answer;
}

double NCPA::mean( double *vec, size_t n ) {
	double running = 0.0;
	for (size_t i = 0; i < n; i++) {
		running += vec[i];
	}
	return running / ((double)n);
}

void NCPA::pol2cart( double r, double theta_rad, double &x, double &y ) {
	x = r * std::cos( theta_rad );
	y = r * std::sin( theta_rad );
}


void cart2pol( double x, double y, double &r, double &theta_rad ) {
	r = std::sqrt( x * x + y * y );
	theta_rad = std::atan2( y, x );
}

void NCPA::read_matrix_from_file( const std::string &filename, double **&contents,
	size_t &ncols, size_t &nrows, std::string delimiters ) {

	std::vector< std::vector< std::string > > textcontents;
	read_text_columns_from_file( filename, textcontents, delimiters );

	ncols = textcontents.size();
	nrows = textcontents[0].size();

	// allocate the memory
	contents = NCPA::allocate_matrix<double>( ncols, nrows );

	for (size_t colind = 0; colind < ncols; colind++) {
		std::vector< std::string > *vecPtr = &(textcontents[colind]);
		for (size_t rowind = 0; rowind < nrows; rowind++) {
			contents[colind][rowind] = std::stod( vecPtr->at(rowind) );
		}
	}
}

void NCPA::read_text_columns_from_file( const std::string &filename,
		std::vector< std::vector< std::string > > &contents, std::string delimiters ) {

	// count the rows
	//size_t nrows = count_rows_arbcol( filename );

	// count the columns
	std::ifstream infile( filename );
	std::string line;
	NCPA::safe_getline( infile, line );
	std::vector< std::string > fields = NCPA::split( line, delimiters );
	size_t ncols = fields.size();

	// allocate
	size_t colind;
	contents = std::vector< std::vector< std::string > >( ncols );
	for (colind = 0; colind < ncols; colind++) {
		contents[ colind ].resize( 1 );
		contents[ colind ][ 0 ] = fields[ colind ];
	}

	// read the rest
	while (infile.good()) {
		NCPA::safe_getline( infile, line );
		fields.clear();
		fields = NCPA::split( line, delimiters );
		if (fields.size() == ncols) {
			for (colind = 0; colind < ncols; colind++) {
				contents[colind].push_back( fields[ colind ] );
			}
		}
	}
	infile.close();
}

// void NCPA::linspace( double a, double b, size_t k, double *&ls ) {
// 	double d = (b - a) / ((double)(k - 1));
// 	for (size_t i = 0; i < k; i++) {
// 		ls[ i ] = a + ((double)i) * d;
// 	}
// }

// void NCPA::logspace( double a, double b, size_t k, double *&ls ) {
// 		double la = std::log10(a);
// 		double lb = std::log10(b);
// 		NCPA::linspace( la, lb, k, ls );
// 		for (size_t i = 0; i < k; i++) {
// 			ls[ i ] = std::pow( 10.0, ls[ i ] );
// 		}
// 	}

std::vector<double> NCPA::random_numbers( size_t N_, double scale ) {
	std::vector<double> randn;
	randn.reserve( N_ );
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<> distribution(0.0, scale);
	for (size_t i = 0; i < N_; i++) {
		double r = distribution( generator );
		randn.push_back( r );
	}
	return randn;
}

// @todo Start the count at floor(log2(v)) and test
size_t NCPA::nextpow2( size_t v ) {
	size_t p = 0;
	double dv = (double)v;
	while (std::exp2(p) < dv) {
		p++;
	}
	return p;
}