#include "BroadbandPropagator.h"

#include "util.h"

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef MAX_MODES
#  define MAX_MODES 4000
#endif

#ifndef Pi
#  define Pi 3.141592653589793
#endif

#define FFTN 32 * 1024

NCPA::BroadbandPropagator::~BroadbandPropagator() {
    // delete[] r_vec;
    // delete[] f_vec;
}

// we have established the relationship between FFTW_BACKWARD and Matlab ifft
//      FFTW_BACKWARD(x,NFFT) = NFFT*ifft_matlab(x, NFFT)
// and  FFTW_FORWARD(x,NFFT)  = fft_matlab(x, NFFT)
// Thus spectrum = FFTW_BACKWARD(pulse,NFFT)*dt
//      pulse    = FFTW_FORWARD(spectrum,NFFT)*df
// where df = f_max/f_step and dt = 1/(NFFT*df)
//
// Added NFFT as argument; normalizes the builtin pulse
// the source spectrum is loaded/computed in this function
//
// int NCPA::BroadbandPropagator::get_source_spectrum(
//     std::complex<double> *dft_vec, std::complex<double> *pulse_vec,
//     std::complex<double> *arg_vec ) {
int NCPA::BroadbandPropagator::get_source_spectrum(
    std::vector<std::complex<double>>& dft_vec,
    std::vector<std::complex<double>>& pulse_vec,
    std::vector<std::complex<double>>& arg_vec ) {
    //
    // Note that this function returns the source spectrum at 0 and positive
    // frequencies!
    //
    int i;
    double dt, fmx, scale;
    std::complex<double> I( 0.0, 1.0 );
    FILE *f;
    fftw_plan p;

    scale = 0.0;   // initialize 'scale'
    fmx = ( (double)NFFT )
        * f_step;  // max frequency; the resulting time interval is dt = 1/fmx
    dt  = 1.0 / fmx;

    if (source_type.compare( "impulse" )
        == 0) {  // use spectrum of a delta function => to get the impulse
                 // response;
        for (i = 0; i < Nfreq; i++) {
            dft_vec[ i ]
                = 1.0 + I * 0.0;  // 'ones' for all positive frequencies
            dft_vec[ i ] = dft_vec[ i ]
                         * exp( I * 2.0 * Pi * f_vec[ i ] * (double)Nfreq
                                / 2.0 );  // add delay - DV 20150720
        }
    } else if (source_type.compare( "pulse1" ) == 0
               || source_type.compare( "pulse2" ) == 0) {
        // use one of the two built-in pulses (provided by Roger Waxler)

        // f_center is initialized with negative value if the option --f_center
        // was not used if that's the case redefine f_center to be = f_max/5;
        if (f_center < 0.0) {
            f_center = f_vec[ Nfreq - 1 ] / 5.0;
        }

        // get the scale for the built-in pulse
        if (( f_center > 0 )
            & ( f_center <= f_vec[ Nfreq - 1 ] / 5.0 + 1.0E-16 )) {
            scale = 1.0 / f_center;
        } else {
            std::ostringstream es;
            es << std::endl
               << "f_center = " << f_center << " Hz is too large." << std::endl
               << "For the built-in pulse f_center should be set smaller than "
                  "f_max/5 = "
               << f_vec[ Nfreq - 1 ] / 5.0;
            throw std::invalid_argument( es.str() );
        }

        if (source_type.compare( "pulse1" ) == 0) {
            // make sure the dispersion file contains what we need to generate
            // a good pulse pulse center frequency f_center should be set <=
            // f_max/5


            // get the source spectrum at 0 and positive frequencies
            // divide by 2 since the pulse has amplitude one and we only want
            // the spectrum at positive freqs after propagation and ifft we'll
            // multiply by 2 to get the propagated pulse amplitude
            for (i = 0; i < Nfreq; i++) {
                dft_vec[ i ] = pulse_spec_fit( scale, f_vec[ i ] )
                             * std::exp( I * 2.0 * Pi * f_vec[ i ] * scale )
                             / 2.0;  // source spectrum
            }
        }  // end pulse choice # 1

        // ----------- do pulse choice # 2 ---------------------
        else if (source_type.compare( "pulse2" ) == 0) {
            double power = 8.0;  // 'power' hardcoded here
            // scale = 1.0/f_center;

            // cout << "scale = " << scale << endl;

            // arg_vec is the FT of initial pulse choice #2
            // dft_vec = model_pulse_fft( power, scale, dt );
            dft_vec = model_pulse_fft( power, scale, dt );
            dft_vec.resize( Nfreq );
            // for (i = 0; i < Nfreq; i++) {
            //     //  multiplication by the dt factor
            //     //  already done in model_pulse_fft()
            //     dft_vec[ i ] = arg_vec[ i ];
            // }
        }  // end pulse choice # 2
    }  // end of builtin pulse choices 1 and 2

    else if (source_type.compare( "spectrum" )
             == 0) {  // use custom source spectrum from a file
        std::cout << "Loading source spectrum from file " << source_file
                  << std::endl;
        // double *freqv;
        // freqv = new double[ Nfreq ];
        std::vector<double> freqv( Nfreq );
        load_source_spectrum( freqv, dft_vec );
        if (fabs( freqv[ 0 ] - f_vec[ 0 ] ) > 1.0e-8) {  // compare frequencies
            std::cerr << "The frequencies from the source spectrum file "
                      << source_file
                      << " do not appear to match the ones from the provided "
                         "dispersion file"
                      << std::endl
                      << " ... aborting." << std::endl;
            exit( 1 );
        }
        // delete[] freqv;
    }

    if (source_type.compare( "waveform" ) == 0) {
        std::vector<double> t, tdp;
        load_source_pulse_td( t, tdp );
        dt = t[ 1 ] - t[ 0 ];
        for (i = 0; i < (int)( t.size() ); i++) {
            pulse_vec[ i ]
                = tdp[ i ];  // automatic conversion to complex<double>
        }
        //
        // perform integral of (pulse(t))*exp(+iwt)*dt via fft
        // to obtain the source spectrum (with NFFT points): 'arg_vec'
        // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt
        // factor afterwards Note that there is no need to multiply by N as in
        // Matlab; note: this FFTW_BACKWARD lacks the 1/N factor that Matlab
        // ifft() has in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
        //
        p = fftw_plan_dft_1d(
            NFFT, reinterpret_cast<fftw_complex *>( pulse_vec.data() ),
            reinterpret_cast<fftw_complex *>( arg_vec.data() ), FFTW_BACKWARD,
            FFTW_ESTIMATE );
        fftw_execute( p );
        fftw_destroy_plan( p );

        // multiply by the dt factor to complete the Fourier integral over time
        // note: it is expected that the energy in the pulse is concentrated
        // well below f_max
        for (i = 0; i < Nfreq; i++) {
            dft_vec[ i ] = arg_vec[ i ] * dt;
        }
    } else {
        for (i = 0; i * f_step < 0.99 * f_vec[ 0 ]; i++)
            arg_vec[ i ]
                = 0.0;  // left zero pad up to f_min present in the spectrum;
        int i0 = i;
        for (i = i0; i < ( i0 + Nfreq ); i++)
            arg_vec[ i ] = dft_vec[ i - i0 ];  // arg_vec contains src spectrum
        for (; i < NFFT; i++)
            arg_vec[ i ]
                = 0.0;  // arg_vec[] zero pad to the right of the src spectrum

        //
        // perform fft to obtain the pulse at the source (time domain) of NFFT
        // points: 'pulse_vec' Note though that since arg_vec contains only the
        // positive freq spectrum we'll have to double the result when we
        // convert to the time domain
        //
        p = fftw_plan_dft_1d(
            NFFT, reinterpret_cast<fftw_complex *>( arg_vec.data() ),
            reinterpret_cast<fftw_complex *>( pulse_vec.data() ), FFTW_FORWARD,
            FFTW_ESTIMATE );
        fftw_execute( p );
        fftw_destroy_plan( p );

        // multiply by f_step to complete the Fourier integral
        for (i = 0; i < NFFT; i++) {
            pulse_vec[ i ] = f_step * pulse_vec[ i ];
        }
    }

    // save initial pulse: pulse_vec
    if (1) {
        f = fopen( "source_waveform_input.dat", "w" );
        if (source_type.compare( "pulse1" ) == 0
            || source_type.compare( "pulse2" )
                   == 0) {  // the builtin pulse choices 1 and 2
            double factor = 2.0;
            // for(i=0;i<(NFFT*f_step*5*scale);i++) {
            for (i = 0; i < NFFT / 2; i++) {
                // fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx,
                // real(pulse_vec[i]), imag(pulse_vec[i]));
                //  save 2*real(pulse_vec[i]); factor of 2 because we only had
                //  the spectrum for positive freqs
                fprintf( f, "%12.6f %15.6e\n", 1.0 * i / fmx,
                         factor * real( pulse_vec[ i ] ) );
                // printf("%12.6f %15.6e\n", 1.0*i/fmx,
                // factor*real(pulse_vec[i]));
            }
        } else {  // if not the builtin pulse
            double factor = 2;
            if (source_type.compare( "waveform" ) == 0)
                factor
                    = 1;  // use custom source pulse (time domain) from a file

            for (i = 0; i < NFFT / 2; i++) {
                // fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx,
                // real(pulse_vec[i]), imag(pulse_vec[i]));
                //  save 2*real(pulse_vec[i]); factor of 2 if we only had
                //  spectrum for positive freqs
                fprintf( f, "%12.6f %15.6e\n", 1.0 * i / fmx,
                         factor * real( pulse_vec[ i ] ) );
            }
        }
        fclose( f );
        // printf("src_flg = %d\n", src_flg);
        std::cout << "--> Initial source waveform saved in file "
                     "'source_waveform_input.dat' "
                  << std::endl
                  << " with format: | Time (s) | Amplitude |" << std::endl;
    }

    // save the source spectrum
    f = fopen( "source_spectrum_input.dat",
               "w" );  // save the source spectrum in this file
    for (i = 0; i < Nfreq; i++) {
        fprintf( f, "%14.9f %15.6e %15.6e\n", f_vec[ i ], real( dft_vec[ i ] ),
                 imag( dft_vec[ i ] ) );
    }
    fclose( f );
    std::cout << "--> The source spectrum "
              << "is saved in file 'source_spectrum_input.dat' " << std::endl
              << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << std::endl;

    return 0;
}  // end of get_source_spectrum()

std::complex<double> NCPA::BroadbandPropagator::pulse_spec_fit( double scale,
                                                                double x ) {
    // double fnorm = 1.0;
    double fnorm
        = 0.780776406404415;  // this is the value ((sqrt(4.0+0.25)-0.5)*0.5);
    std::complex<double> answer, I;
    I = std::complex<double>( 0.0, 1.0 );

    /* scale down the frequency to get infrasound */

    x = scale * fnorm * x;

    /* fit to fourier transform of signal at 10 m */
    answer = x * exp( ( -x - x * x ) / 2 ) * exp( I * ( -0.5 * Pi ) ) * scale
           / 33.4336;
    /* calibrate (times 10 divided by 2 for pressure doubling) and time shift
     * by 4*scale sec */
    /* multiply by 4 pi to get a delta func source strength */
    answer = 4.0 * Pi * 5.0 * answer * exp( I * 2.0 * Pi * x * ( 0.01 ) );

    return answer;
}

// void NCPA::BroadbandPropagator::model_pulse_fft(
//     double power, double scale, double dt, std::complex<double> *target_vec
//     ) {
std::vector<std::complex<double>> NCPA::BroadbandPropagator::model_pulse_fft(
    double power, double scale, double dt ) {
    int npoints = NFFT;
    std::vector<std::complex<double>> pulse( npoints );
    std::vector<std::complex<double>> spectrum( npoints );

    for (int i = 0; i < npoints; i++) {
        // the pulse timeseries
        pulse[ i ].real( model_pulse_shape( power, scale, i * dt ) );
    }

    // Perform Fourier transform to obtain the spectrum of the initial pulse
    fftw_plan p = fftw_plan_dft_1d(
        npoints, reinterpret_cast<fftw_complex *>( pulse.data() ),
        reinterpret_cast<fftw_complex *>( spectrum.data() ), FFTW_BACKWARD,
        FFTW_ESTIMATE );

    fftw_execute( p );
    fftw_destroy_plan( p );

    // multiply by dt to complete the Fourier integral
    for (int i = 0; i < npoints; i++) {
        spectrum[ i ] *= dt;
    }
    return spectrum;
}

/* Power law pulse model with exponential decay. */
double NCPA::BroadbandPropagator::model_pulse_shape( double power,
                                                     double scale, double x ) {
    double answer;

    x      = 2.094 * ( 1.0 + power ) * x / scale;
    answer = pow( x, power ) * ( 1 - x / ( 1.0 + power ) ) * exp( -x )
           / 1.387566e+03;

    return answer;
}

// int NCPA::BroadbandPropagator::load_source_spectrum(
//     double *freqv, std::complex<double> *dft_vec ) {
int NCPA::BroadbandPropagator::load_source_spectrum(
    std::vector<double>& freqv, std::vector<std::complex<double>>& dft_vec ) {
    // load source spectrum
    int i;
    double d1, d2, d3;
    std::complex<double> I( 0.0, 1.0 );
    std::ifstream indata;

    indata.open( source_file.c_str() );
    if (!indata) {
        std::cerr << "file " << source_file << " could not be opened.\n";
        exit( 1 );
    }

    i = 0;
    indata >> d1 >> d2 >> d3;  // read first line in file
    while (!indata.eof()) {
        freqv[ i ]   = d1;
        dft_vec[ i ] = d2 + I * d3;
        indata >> d1 >> d2 >> d3;
        if (i > ( NFFT - 1 )) {
            std::cerr << "file " << source_file
                      << "appears to have more points than maximum allowed of "
                      << std::endl
                      << "NFFT = " << NFFT << std::endl
                      << " ... aborting." << std::endl;
            exit( 1 );
        }
        i++;
    }
    indata.close();

    return 0;
}

int NCPA::BroadbandPropagator::load_source_pulse_td(
    std::vector<double>& t, std::vector<double>& tdp ) {
    // load time domain source pulse: time (s) | amplitude
    std::cout << "Loading time domain source pulse from file " << source_file
              << std::endl;

    double d1, d2;  //, d3;
    std::ifstream indata;

    indata.open( source_file.c_str() );
    if (!indata) {
        std::cerr << "file " << source_file << " could not be opened."
                  << std::endl;
        exit( 1 );
    }

    indata >> d1 >> d2;  // read first line in file
    while (!indata.eof()) {
        t.push_back( d1 );
        tdp.push_back( d2 );
        indata >> d1 >> d2;
    }
    indata.close();

    return 0;
}

// void NCPA::BroadbandPropagator::fft_pulse_prop(
//     double t0, double range, std::complex<double> *dft_vec,
//     std::complex<double> *pulse_vec ) {
void NCPA::BroadbandPropagator::fft_pulse_prop(
    double t0, double range, std::vector<std::complex<double>>& dft_vec,
    std::vector<std::complex<double>>& pulse_vec ) {
    int i, i0;  //,j,smooth_space;
    // double sqrt_rho_ratio = sqrt(rho_zrcv/rho_zsrc);
    std::complex<double> t_phase;
    std::complex<double> I( 0.0, 1.0 );
    // std::complex<double> expov8pir = I*exp(-I*Pi*0.25)/sqrt(8.0*Pi*range);
    // double df = f_vec[ 1 ] - f_vec[ 0 ];
    double df = f_vec[ 2 ] - f_vec[ 1 ];  // because the DC is always included
    fftw_plan p;
    std::cout << "calculated df = " << df << std::endl;
    if (NFFT < Nfreq) {
        throw std::invalid_argument(
            "fft too short (i.e. NFFT < Nfreq), exiting." );
    }

    // arg_vec = new std::complex<double>[ NFFT ];
    std::vector<std::complex<double>> arg_vec( NFFT, 0.0 );

    // left zero pad
    for (i = 0; i * df < 0.99 * f_vec[ 0 ]; i++)
        arg_vec[ i ] = 0.0;  // left zero pad

    i0 = i;
    for (i = 0; i < Nfreq; i++) {
        // calculate phase shift
        // corresponding to reduced time t0; note f_vec[i+1]
        t_phase = exp( -I * 2.0 * Pi * f_vec[ i ] * t0 );

        // the reduced pressure is defined as:
        // (actual pressure)/sqrt(rho(z)): p_reduced(r,z) = p(r,z)/sqrt(rho(z))

        // Note on mode scaling in this code vs. the modes in Computational Oc.
        // Acoust. 1994: V_book =  sqrt(rho)*V_in_this_code; Thus the formula
        // for the Fourier component of pressure using the modes in this code
        // is given in DV Modess notes eq. 25 pg. 3 and contains the factor
        // sqrt_rho_ratio as opposed to the factor 1/rho(z_s) in the book
        // eq. 5.14 pg 274. will multiply by df at the end; note f_vec[i+1]
        arg_vec[ i0 + i ] = t_phase * dft_vec[ i ]
                          * transfer_function[ i ];  // note sqrt_rho_ratio
    }

    for (; i < NFFT; i++) arg_vec[ i ] = 0.0;  // right zero pad

    //
    // perform fft of arg_vec to obtain the propagated time domain pulse
    //
    p = fftw_plan_dft_1d( NFFT,
                          reinterpret_cast<fftw_complex *>( arg_vec.data() ),
                          reinterpret_cast<fftw_complex *>( pulse_vec.data() ),
                          FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute( p );
    fftw_destroy_plan( p );

    // multiply by df to complete the Fourier integral
    for (i = 0; i < NFFT; i++) {
        pulse_vec[ i ] *= df;
    }

    // delete[] arg_vec;
}  // end of debug version of 'fft_pulse_prop'

double NCPA::BroadbandPropagator::half_hann( int begin, int end, int i ) {
    double answer;
    if (i < begin)
        answer = 1.0;
    else if (( begin <= i ) && ( i <= end )) {
        answer = 0.5
               * ( std::cos( Pi * (( i - begin )) / ( (int)( end - begin ) ) )
                   + 1.0 );
    } else
        answer = 0.0;
    return answer;
}
