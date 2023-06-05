#ifndef NCPA_FFT_H_INCLUDED
#define NCPA_FFT_H_INCLUDED

#include <complex>
#include "fftw3.h"

namespace NCPA {

	// FFTW interface
	void fft( size_t N, double *in, std::complex<double> *&out );
	void fft( size_t N, std::complex<double> *in, std::complex<double> *&out );
	void ifft( size_t N, double *in, std::complex<double> *&out );
	void ifft( size_t N, std::complex<double> *in, std::complex<double> *&out );
	void do_transform( size_t N, int direction, fftw_complex *in,
		fftw_complex *&out );
}



#endif
