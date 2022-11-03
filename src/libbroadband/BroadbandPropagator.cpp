#include "BroadbandPropagator.h"
#include <complex>
#include <fftw3.h>
#include <string>
#include <stdexcept>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include "util.h"

#ifndef MAX_MODES
#define MAX_MODES 4000
#endif

#ifndef Pi
#define Pi 3.141592653589793
#endif

#define FFTN 32*1024

// NCPA::BroadbandPropagator::BroadbandPropagator( NCPA::ParameterSet *param ) {

// 	// read in and store parameters
// 	NFFT 					= param->getInteger( "nfft" );
// 	waveform_out_file 		= param->getString( "output_waveform_file" );
// 	dispersion_input_file 	= param->getString( "input_dispersion_file" );
// 	source_type 			= param->getString( "source" );
// 	source_file 			= param->getString( "source_file" );
// 	f_center 				= param->getFloat( "f_center" );
// 	max_cel 				= param->getFloat( "max_celerity" );
// 	single_receiver 		= param->getString("receiver").compare("single") == 0;

// 	if (single_receiver) {
// 		Nr = 1;
// 		r_vec = new double[ 1 ];
// 		r_vec[ 0 ] = param->getFloat( "range_km" ) * 1000.0;
// 	} else {
// 		double rmin, rmax, rstep;
// 		rmin = param->getFloat( "start_range_km" );
// 		rmax = param->getFloat( "end_range_km" );
// 		rstep = param->getFloat( "range_step_km" );
// 		if (rmax <= rmin) {
// 			throw std::runtime_error( "end_range_km must be greater than start_range_km" );
// 		}
// 		Nr = (int)floor( (rmax-rmin) / rstep ) + 1;
// 		r_vec = new double[ Nr ];
// 		for (int i = 0; i < Nr; i++) {
// 			r_vec[ i ] = (rmin + ((double)i) * rstep) * 1000.0;  // expected in meters
// 		}
// 	}

// 	Nfreq = NCPA::count_rows_arbcol( dispersion_input_file ) + 1;    // count DC frequency
// 	f_vec = new double[ Nfreq ];
// 	mode_count = new int[ Nfreq ];
// 	re_k       = dmatrix(Nfreq, MAX_MODES);
// 	im_k       = dmatrix(Nfreq, MAX_MODES);
// 	mode_S     = dmatrix(Nfreq, MAX_MODES);
// 	mode_R     = dmatrix(Nfreq, MAX_MODES);
// 	read_dispersion_file();

// 	f_step = f_vec[ 1 ] - f_vec[ 0 ];
// }


NCPA::BroadbandPropagator::~BroadbandPropagator() {
	delete [] r_vec;
	delete [] f_vec;
	// delete [] mode_count;
	// free_dmatrix(re_k,   Nfreq, MAX_MODES);
	// free_dmatrix(im_k,   Nfreq, MAX_MODES);								
	// free_dmatrix(mode_S, Nfreq, MAX_MODES);
	// free_dmatrix(mode_R, Nfreq, MAX_MODES);
}

// int NCPA::BroadbandPropagator::read_dispersion_file() {

//   // this function reads the file created with write_dispersion()
//   // reads freq #modes rho(@z=src), rho(@z=rcv), k_horiz, modes(@z=src), modes(@z=receiverheight)
//   int m,n;
//   double x1,x2,y1,y2;
//   FILE *fp = fopen(dispersion_input_file.c_str(),"r");
//   if(fp==NULL){
//       std::ostringstream es;
//       es << "file << " << dispersion_input_file.c_str() << " could not be opened.\n";
//       throw std::invalid_argument(es.str());
//   }
//   //read_header(f); // use it if file has header
//   printf("--> Opening dispersion file %s ...if the program chokes here the file may be corrupt.\n",
//                              dispersion_input_file.c_str());
//   // read line by line;
//   f_vec[0] = 0; // the DC frequency
//   mode_count[0] = 0;
//   n=1; // start from 1 because f_vec[0] = 0
//   //n = 0;
//   while(fscanf(fp,"%lf",&f_vec[n])==1){
//       fscanf(fp,"%d",&mode_count[n]);
//       fscanf(fp,"%lf",&rho_zsrc);
//       fscanf(fp,"%lf",&rho_zrcv);
//       for(m=0;m<mode_count[n];m++){
//           if(fscanf(fp,"%lf %lf %lf %lf",&x1,&x2,&y1,&y2)!=4){
//               fprintf(stderr,"%s: read_dispersion_data format error\n",dispersion_input_file.c_str());
//               exit(1);
//           }
//           re_k[n][m]  =x1;
//           im_k[n][m]  =x2;
//           mode_S[n][m]=y1;
//           mode_R[n][m]=y2;
//       }
//       ++n;
//   }
//   fclose(fp);
//   printf("--> Found modes at %d positive frequencies in file %s\n",n, dispersion_input_file.c_str());

//   return 0; 
// }

// int NCPA::BroadbandPropagator::calculate() {
// 	int i,n;
//   double rr, tskip, fmx, t0;	
//   std::complex<double> cup, *dft_vec, *pulse_vec, *arg_vec;
//   //complex<double> I = complex<double> (0.0, 1.0);
//   FILE *f;
  
//   // if the option --nfft was not passed to the main program then the 
//   // NFFT default was -1. Otherwise it has the requested positive value. 
//   // Here we make sure that whatever the current value of NFFT is 
//   // it's not less than 4*n_freqs
//   if (NFFT < Nfreq) {
//     NFFT = 4 * Nfreq;
//     printf("Minimum NFFT is set to NFFT = %d\n", NFFT);
//   }
  
  
//   // if zero attenuation requested set the imaginary part of the wavenumber = 0
//   if (zero_attn_flag) {
//     for (n=0; n<Nfreq; n++) {
//       for(int m=0;m<mode_count[n];m++){
//           im_k[n][m]  = 0.0;
//       }
//     }
//   }

//   dft_vec   = new std::complex<double> [Nfreq];
//   modal_sum = new std::complex<double> [Nfreq];
//   pulse_vec = new std::complex<double> [NFFT];
//   arg_vec   = new std::complex<double> [NFFT]; 
//   std::memset( dft_vec,   0, Nfreq * sizeof( std::complex< double > ) );
//   std::memset( modal_sum, 0, Nfreq * sizeof( std::complex< double > ) );
//   std::memset( pulse_vec, 0, NFFT  * sizeof( std::complex< double > ) );
//   std::memset( arg_vec,   0, NFFT  * sizeof( std::complex< double > ) );

//   fmx = ((double)NFFT)*f_step; // max frequency
								      
//   get_source_spectrum( dft_vec, pulse_vec, arg_vec );
	
// 	// if (Nr == 1) { // propagation to one receiver at distance RR from source
// 	//     std::cout << "--> Doing pulse propagation source-to-receiver at one range: " 
// 	//          << r_vec[0]/1000.0 << " km" << std::endl;
	         
// 	//     t0 = r_vec[0] /max_cel;
//  //      std::memset( modal_sum,  0, Nfreq * sizeof( std::complex< double > ) );
//  //      compute_modal_sum( r_vec[ 0 ] );
// 	//     //
// 	//     // fft propagation: note that here 'pulse_vec' is overwritten with the propagated pulse
// 	//     //							
//  //      fft_pulse_prop( t0, r_vec[0], dft_vec, pulse_vec);						      
      
//  //      // save propagated pulse to file
// 	//   // DV 20150514 - factor of two necessary because we only used 
// 	//   // the positive freq. spectrum in the fft
//  //      // DV 20170810 - parameter factor to make it easy to agree with other codes 
//  //      // (e.g. Roger Waxler's modal code)
//  //      double factor = 2.0;													
//  //      f = fopen(waveform_out_file.c_str(),"w");
//  //      for(i=0; i<NFFT; i++) {
//  //          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx+t0, factor*real(pulse_vec[i]));
//  //      }
//  //      fclose(f);

//  //      std::cout << "--> Propagated pulse saved in file: '" << waveform_out_file << "'" << std::endl
//  //           << "' with format: | Time (s) | Re(pulse) |" << std::endl;  

// 	// }
// 	// else { // propagation to several receivers
// 		// @todo combine these, only difference is the use of multiple ranges?

// 	  std::cout << "--> Propagating pulse from source-to-receivers on grid ..." << std::endl;					      
//       printf("----------------------------------------------\n");
//       printf("max_celerity     t0             R\n");
//       printf("    m/s          sec            km\n");
//       printf("----------------------------------------------\n");

//       f=fopen(waveform_out_file.c_str(),"w");
//       tskip = 0.0;
//       //double DR = r_vec[1] - r_vec[0];
//       //for(n=0; n<=(int)(std::floor((r_vec[Nr-1]-r_vec[0])/DR)); n++) {
//       for (n=0; n < Nr; n++) {
//           //rr= R_start + DR*n;
//       	  rr = r_vec[ n ];
//           t0=tskip+rr/max_cel;
//           printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

//           std::memset( modal_sum,  0, Nfreq * sizeof( std::complex< double > ) );
//           compute_modal_sum( rr );

//           // the call with NFFT as argument			            
//           fft_pulse_prop( t0, rr, dft_vec, pulse_vec );

//           // DV 20170810 - parameter 'factor' to make it easy to agree with other codes 
//           // (e.g. Roger Waxler's modal code)
//           double factor = 2.0;								              
//           for(i=0;i<NFFT;i++){
//               fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx+t0, factor*real(pulse_vec[i]));
//           }
//           fprintf(f,"\n");
//       }
//       fclose(f);
//       printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
//       printf("Time array length = %d; delta_t = %g s\n", NFFT, 1.0/fmx);
//       printf("Propagation results saved in file: %s\n", waveform_out_file.c_str());
//       printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
//   // }

//   delete [] pulse_vec;
//   delete [] arg_vec;
//   delete [] dft_vec;
//   delete [] modal_sum;

//   return 0;
// }


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
int NCPA::BroadbandPropagator::get_source_spectrum( std::complex<double> *dft_vec, 
	std::complex<double> *pulse_vec, std::complex<double> *arg_vec ) 
{
  //
  // Note that this function returns the source spectrum at 0 and positive frequencies!
  //
  int i;
  double dt, fmx, scale;
  std::complex<double> I = std::complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  scale = 0.0; // initialize 'scale'
  fmx = ((double)NFFT)*f_step; // max frequency; the resulting time interval is dt = 1/fmx
  dt  = 1.0/fmx;
  
  if (source_type.compare("impulse") == 0) { // use spectrum of a delta function => to get the impulse response;
      for(i=0; i<Nfreq; i++) {
          dft_vec[i] = 1.0 + I*0.0; // 'ones' for all positive frequencies
          dft_vec[i] = dft_vec[i]*exp(I*2.0*Pi*f_vec[i]*(double)Nfreq/2.0); // add delay - DV 20150720
      }
  } else if ( source_type.compare("pulse1") == 0 || source_type.compare("pulse2") == 0 ) { 
  //use one of the two built-in pulses (provided by Roger Waxler)

  	  // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be = f_max/5;
      if (f_center < 0.0) {
          f_center = f_vec[Nfreq-1] / 5.0;
      }

      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[Nfreq-1]/5.0+1.0E-16)) {	
          scale = 1.0/f_center;
      } else {
          std::ostringstream es;
          es << std::endl << "f_center = " << f_center << " Hz is too large." << std::endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[Nfreq-1]/5.0;
          throw std::invalid_argument(es.str());
      }
          
      if (source_type.compare("pulse1") == 0) {  
          // make sure the dispersion file contains what we need to generate a good pulse
          // pulse center frequency f_center should be set <= f_max/5
          
          
         
          // get the source spectrum at 0 and positive frequencies
          // divide by 2 since the pulse has amplitude one and we only want the spectrum at positive freqs
          // after propagation and ifft we'll multiply by 2 to get the propagated pulse amplitude
          for(i=0; i<Nfreq; i++) {
              dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale)/2.0; // source spectrum    
          }
      } // end pulse choice # 1

      // ----------- do pulse choice # 2 ---------------------
      else if (source_type.compare("pulse2") == 0) {
          double power = 8.0; // 'power' hardcoded here
          //scale = 1.0/f_center;
        
          //cout << "scale = " << scale << endl;

          model_pulse_fft(power, scale, dt, arg_vec); // arg_vec is the FT of initial pulse choice #2
          for (i=0; i<Nfreq; i++) {
              dft_vec[i] = arg_vec[i]; //  multiplication by the dt factor already done in model_pulse_fft()
          }
      } // end pulse choice # 2   
  } // end of builtin pulse choices 1 and 2

  else if (source_type.compare("spectrum") == 0) { //use custom source spectrum from a file
      std::cout << "Loading source spectrum from file " << source_file << std::endl;
      double *freqv;
      freqv = new double [Nfreq];
      load_source_spectrum(freqv, dft_vec);
      if (fabs(freqv[0]-f_vec[0])>1.0e-8) { // compare frequencies
          std::cerr << "The frequencies from the source spectrum file " << source_file
               << " do not appear to match the ones from the provided dispersion file" 
               << std::endl << " ... aborting." << std::endl;
               exit(1);
      }
      delete [] freqv;
  }

  if (source_type.compare("waveform") == 0) {
      std::vector<double> t, tdp;
      load_source_pulse_td(t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)(t.size()); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            					reinterpret_cast<fftw_complex*> (arg_vec), \
                                FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<Nfreq; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
	} else {
      for(i=0;i*f_step<0.99*f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum;
      int i0 = i;
      for(i=i0; i<(i0+Nfreq); i++) arg_vec[i] = dft_vec[i-i0]; // arg_vec contains src spectrum
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      // Note though that since arg_vec contains only the positive freq spectrum
      // we'll have to double the result when we convert to the time domain
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        				reinterpret_cast<fftw_complex*> (pulse_vec), \
			        							FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);

      // multiply by f_step to complete the Fourier integral
      for(i=0;i<NFFT;i++) {
        pulse_vec[i]=f_step*pulse_vec[i];
      }   
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform_input.dat","w");
      if ( source_type.compare("pulse1") == 0 || source_type.compare("pulse2") == 0 ) { // the builtin pulse choices 1 and 2
        double factor = 2.0;
        //for(i=0;i<(NFFT*f_step*5*scale);i++) {
        for(i=0;i<NFFT/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had the spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, factor*real(pulse_vec[i]));
            //printf("%12.6f %15.6e\n", 1.0*i/fmx, factor*real(pulse_vec[i]));
        }
      } else { // if not the builtin pulse
        double factor = 2;
        if (source_type.compare("waveform") == 0) 
        	factor = 1;  // use custom source pulse (time domain) from a file
        
        for(i=0;i<NFFT/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, factor*real(pulse_vec[i]));
        }    
      }
      fclose(f);
      //printf("src_flg = %d\n", src_flg);
      std::cout << "--> Initial source waveform saved in file 'source_waveform_input.dat' "
      	 << std::endl << " with format: | Time (s) | Amplitude |" << std::endl;
  }

  // save the source spectrum
  f=fopen("source_spectrum_input.dat","w"); // save the source spectrum in this file
  for(i=0; i<Nfreq; i++) {
    fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
  }
  fclose(f);
  std::cout << "--> The source spectrum "
       << "is saved in file 'source_spectrum_input.dat' " << std::endl
       << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << std::endl; 

  return 0;
} // end of get_source_spectrum()


std::complex<double> NCPA::BroadbandPropagator::pulse_spec_fit(double scale,double x){

  //double fnorm = 1.0;
  double fnorm = 0.780776406404415; // this is the value ((sqrt(4.0+0.25)-0.5)*0.5);
  std::complex<double> answer, I;
  I = std::complex<double> (0.0, 1.0);

/* scale down the frequency to get infrasound */

  x=scale*fnorm*x;

  /* fit to fourier transform of signal at 10 m */
  answer=x*exp((-x-x*x)/2)*exp(I*(-0.5*Pi))*scale/33.4336;
  /* calibrate (times 10 divided by 2 for pressure doubling) and time shift by 4*scale sec */
  /* multiply by 4 pi to get a delta func source strength */
  answer=4.0*Pi*5.0*answer*exp(I*2.0*Pi*x*(0.01));

  return answer;
}

void NCPA::BroadbandPropagator::model_pulse_fft(double power,double scale, double dt, 
	std::complex<double> *dft_vec) {
  int i;
  std::complex<double> *p_vec;
  fftw_plan p;

  p_vec = new std::complex<double>( NFFT );

  for(i=0;i<NFFT;i++){
    p_vec[i]=model_pulse_shape(power,scale,i*dt); // the pulse timeseries
  }

  // Perform Fourier transform to obtain the spectrum of the initial pulse
  p=fftw_plan_dft_1d(NFFT, reinterpret_cast<fftw_complex*>(p_vec),
                           reinterpret_cast<fftw_complex*>(dft_vec),
                           FFTW_BACKWARD,FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);

  // multiply by dt to complete the Fourier integral
  for(i=0;i<NFFT;i++){
    dft_vec[i]=dt*dft_vec[i];
  }

  delete [] p_vec;
}

/* Power law pulse model with exponential decay. */
double NCPA::BroadbandPropagator::model_pulse_shape(double power,double scale,double x) {
  double answer;

  x=2.094*(1.0+power)*x/scale;
  answer=pow(x,power)*(1-x/(1.0+power))*exp(-x)/1.387566e+03;

  return answer;
}

int NCPA::BroadbandPropagator::load_source_spectrum(double *freqv, std::complex<double> *dft_vec) {
  // load source spectrum
  int i;
  double d1, d2, d3;
  std::complex<double> I (0.0, 1.0);
  std::ifstream indata;

  indata.open(source_file.c_str());	
  if (!indata) {
      std::cerr << "file " << source_file << " could not be opened.\n";
      exit(1);
  }

  i = 0;
  indata >> d1 >> d2 >> d3; // read first line in file
  while (!indata.eof()) {
      freqv[i] = d1;
      dft_vec[i] = d2 + I*d3;
      indata >> d1 >> d2 >> d3;
      if (i>(NFFT-1)) {
          std::cerr << "file " << source_file 
          		<< "appears to have more points than maximum allowed of " << std::endl
               << "NFFT = " << NFFT << std::endl
               << " ... aborting."  << std::endl;
          exit(1);
      }
      i++;
  }
  indata.close();

  return 0;
}	  



int NCPA::BroadbandPropagator::load_source_pulse_td(std::vector<double> &t, 
	std::vector<double> &tdp ) {

  // load time domain source pulse: time (s) | amplitude
  std::cout << "Loading time domain source pulse from file " << source_file << std::endl;

  double d1, d2; //, d3;
  std::ifstream indata;

  indata.open(source_file.c_str());	
  if (!indata) {
      std::cerr << "file " << source_file << " could not be opened." << std::endl;
      exit(1);
  }

  indata >> d1 >> d2; // read first line in file
  while (!indata.eof()) {
      t.push_back(d1);
      tdp.push_back(d2);
      indata >> d1 >> d2;
  }
  indata.close();

  return 0;
}

// void NCPA::BroadbandPropagator::compute_modal_sum( double range ) {
//   int i, smooth_space, j;
//   std::complex<double> I (0.0, 1.0);
//   std::complex< double > k_H, cup, pot, mode_prod;
//   double sqrt_rho_ratio = sqrt(rho_zrcv/rho_zsrc);
//   std::complex<double> expov8pir = I*exp(-I*Pi*0.25)/sqrt(8.0*Pi*range);

//   for(i=0;i<Nfreq;i++){
//       cup=0.0;
//       for(j=0;j<mode_count[i];j++){
//           k_H = re_k[i][j]+I*im_k[i][j];
//           pot = exp(I*range*k_H);
//           mode_prod = mode_S[i][j]*mode_R[i][j];
//           pot = mode_prod*pot;
//           pot = pot/sqrt(k_H);
//           cup = cup+pot; // modal sum: sum( exp(ikr)/sqrt(k)*Vr*Vs )
//       }
//       // Note on mode scaling in this code vs. the modes in Computational Oc. Acoust. 1994: 
//       // V_book =  sqrt(rho)*V_in_this_code; 
//       // Thus the formula for the Fourier component of pressure using the modes in this code
//       // is given in DV Modess notes eq. 25 pg. 3 and contains the factor sqrt_rho_ratio
//       // as opposed to the factor 1/rho(z_s) in the book eq. 5.14 pg 274.
//       modal_sum[i] = expov8pir*cup*sqrt_rho_ratio; // note sqrt_rho_ratio
//       //arg_vec[i0+i]=sqrt(I/(8.0*Pi*range))*cup;

//   }

//   smooth_space=(int)floor(0.1*Nfreq); // smoothly zero out on right; as in RW (July 2012)

//   for(i=Nfreq-smooth_space;i<Nfreq;i++){
//       modal_sum[i]=modal_sum[i]*half_hann(Nfreq-smooth_space,Nfreq-1,i); // changed df to 1 as df doesn't make sense here
//   }
// }

void NCPA::BroadbandPropagator::fft_pulse_prop( double t0, double range, 
	std::complex<double> *dft_vec, std::complex<double> *pulse_vec ) {

	int i,i0; //,j,smooth_space;
  // double sqrt_rho_ratio = sqrt(rho_zrcv/rho_zsrc);
  std::complex<double> cup,pot,k_H,mode_prod,t_phase,*arg_vec;
  std::complex<double> I (0.0, 1.0);
  // std::complex<double> expov8pir = I*exp(-I*Pi*0.25)/sqrt(8.0*Pi*range);
  double df = f_vec[1] - f_vec[0];
  fftw_plan p;

  if(NFFT < Nfreq){
      throw std::invalid_argument("fft too short (i.e. NFFT < Nfreq), exiting.");
  }
  
  arg_vec = new std::complex<double> [NFFT];

  for(i=0;i*df<0.99*f_vec[0];i++) 
  	arg_vec[i]=0.0; // left zero pad

  i0 = i;
  for(i=0;i<Nfreq;i++){

      // calculate phase shift
      t_phase=exp(-I*2.0*Pi*f_vec[i]*t0); // corresponding to reduced time t0; note f_vec[i+1]
      // cup=0.0;
      // for(j=0;j<mode_count[i];j++){
      //     k_H = re_k[i][j]+I*im_k[i][j];
      //     pot = exp(I*range*k_H);
      //     mode_prod = mode_S[i][j]*mode_R[i][j];
      //     pot = mode_prod*pot;
      //     pot = pot/sqrt(k_H);
      //     cup = cup+pot; // modal sum: sum( exp(ikr)/sqrt(k)*Vr*Vs )
      // }
      // cup=cup*t_phase;
      
      // up to a factor ((delayed Fourier pressure component)*source_spectrum*df) 
      // (see eq. 5.14 pg 274 and eq. 8.9 pg 480 in Comp. Oc. Acoust. 1994 ed.)
      //cup=cup*dft_vec[i]*df;
      // cup=cup*dft_vec[i]; // will multiply by df at the end; note f_vec[i+1]

      // the reduced pressure is defined as: 
      // (actual pressure)/sqrt(rho(z)): p_reduced(r,z) = p(r,z)/sqrt(rho(z))

      // Note on mode scaling in this code vs. the modes in Computational Oc. Acoust. 1994: 
      // V_book =  sqrt(rho)*V_in_this_code; 
      // Thus the formula for the Fourier component of pressure using the modes in this code
      // is given in DV Modess notes eq. 25 pg. 3 and contains the factor sqrt_rho_ratio
      // as opposed to the factor 1/rho(z_s) in the book eq. 5.14 pg 274.
      // will multiply by df at the end; note f_vec[i+1]
      arg_vec[i0+i] = t_phase*dft_vec[i]*transfer_function[i]; // note sqrt_rho_ratio
      //arg_vec[i0+i]=sqrt(I/(8.0*Pi*range))*cup;

  }

  // smooth_space=(int)floor(0.1*Nfreq); // smoothly zero out on right; as in RW (July 2012)

  // for(i=Nfreq-smooth_space;i<Nfreq;i++){
  //     //arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-df,i); // old code : has df in it
  //     arg_vec[i]=arg_vec[i]*half_hann(Nfreq-smooth_space,Nfreq-1,i); // changed df to 1 as df doesn't make sense here
  // }
  for(;i<NFFT;i++) arg_vec[i] = 0.0; // right zero pad
  
  // 
  //perform fft of arg_vec to obtain the propagated time domain pulse
  //
  p=fftw_plan_dft_1d(NFFT, reinterpret_cast<fftw_complex*>(arg_vec), \
                           reinterpret_cast<fftw_complex*>(pulse_vec), \
											FFTW_FORWARD,FFTW_ESTIMATE); 												 
  fftw_execute(p); 
  fftw_destroy_plan(p);
  
  // multiply by df to complete the Fourier integral 
  for(i=0;i<NFFT;i++){
    pulse_vec[i] = df*pulse_vec[i];
  }

  delete [] arg_vec;
}  // end of debug version of 'fft_pulse_prop'


double NCPA::BroadbandPropagator::half_hann(int begin,int end,int i) {
  double answer;
  if(i<begin) answer=1.0;
  else if((begin<=i) && (i<=end)){
      answer=0.5*(std::cos(Pi*((i-begin))/((int)(end-begin)))+1.0);
  }
  else answer=0.0;
  return answer;
}
