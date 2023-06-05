#include "ncpa_fft.h"
#include "fftw3.h"
#include "util.h"
#include <complex>
#include <cstdlib>


void NCPA::fft( size_t N, double *in, std::complex<double> *&out ) {

	fftw_complex *integrand =
		(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	fftw_complex *result = reinterpret_cast<fftw_complex *>( out );
	size_t i;

	// convert from double to fftw_complex
	for (i = 0; i < N; i++) {
		integrand[i][0] = in[ i ];
	}

	NCPA::do_transform( N, FFTW_FORWARD, integrand, result );

	fftw_free(integrand);
}


void NCPA::fft( size_t N, std::complex<double> *in,
		std::complex<double> *&out ) {

	fftw_complex *integrand = reinterpret_cast<fftw_complex *>(in);
	fftw_complex *result = reinterpret_cast<fftw_complex *>(out);

	NCPA::do_transform( N, FFTW_FORWARD, integrand, result );

}


void NCPA::ifft( size_t N, double *in, std::complex<double> *&out ) {
	fftw_complex *integrand =
		(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	fftw_complex *result = reinterpret_cast<fftw_complex *>(out);
	size_t i;

	// convert from double to fftw_complex
	for (i = 0; i < N; i++) {
		integrand[i][0] = in[ i ];
	}

	NCPA::do_transform( N, FFTW_BACKWARD, integrand, result );

	double factor = 1.0 / ((double)N);
	NCPA::vector_scale( N, out, factor, out );

	fftw_free(integrand);
}


void NCPA::ifft( size_t N, std::complex<double> *in,
		std::complex<double> *&out ) {

	fftw_complex *integrand = reinterpret_cast<fftw_complex *>(in);
	fftw_complex *result = reinterpret_cast<fftw_complex *>(out);

	NCPA::do_transform( N, FFTW_BACKWARD, integrand, result );
	double factor = 1.0 / ((double)N);
	NCPA::vector_scale( N, out, factor, out );
}



void NCPA::do_transform( size_t N, int direction,
	fftw_complex *in, fftw_complex *&out ) {

	fftw_plan plan = fftw_plan_dft_1d(N, in, out, direction, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}