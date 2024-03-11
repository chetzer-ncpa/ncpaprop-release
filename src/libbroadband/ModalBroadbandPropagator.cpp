#include "ModalBroadbandPropagator.h"
#include <complex>
#include <fftw3.h>
#include <string>
#include <stdexcept>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <algorithm>
#include "util.h"
#include "parameterset.h"

#ifndef MAX_MODES
#define MAX_MODES 4000
#endif

#ifndef Pi
#define Pi 3.141592653589793
#endif

#define FFTN 32*1024


NCPA::ModalBroadbandPropagator::ModalBroadbandPropagator( NCPA::ParameterSet *param ) {

	// read in and store parameters
	NFFT 					= param->getInteger( "nfft" );
	waveform_out_file 		 = param->getString( "output_waveform_file" );
	dispersion_input_file 	= param->getString( "input_dispersion_file" );
	source_type 			= param->getString( "source" );
	source_file 			= param->getString( "source_file" );
	f_center 				= param->getFloat( "f_center" );
	max_cel 				= param->getFloat( "max_celerity" );
	single_receiver 		= param->getString("receiver").compare("single") == 0;

	if (single_receiver) {
		Nr = 1;
		r_vec = new double[ 1 ];
		r_vec[ 0 ] = param->getFloat( "range_km" ) * 1000.0;
	} else {
		double rmin, rmax, rstep;
		rmin = param->getFloat( "start_range_km" );
		rmax = param->getFloat( "end_range_km" );
		rstep = param->getFloat( "range_step_km" );
		if (rmax <= rmin) {
			throw std::runtime_error( "end_range_km must be greater than start_range_km" );
		}
		Nr = (int)floor( (rmax-rmin) / rstep ) + 1;
		r_vec = new double[ Nr ];
		for (int i = 0; i < Nr; i++) {
			r_vec[ i ] = (rmin + ((double)i) * rstep) * 1000.0;  // expected in meters
		}
	}

	Nfreq = NCPA::count_rows_arbcol( dispersion_input_file ) + 1;    // count DC frequency
	f_vec = new double[ Nfreq ];
	mode_count = new int[ Nfreq ];
	re_k       = dmatrix(Nfreq, MAX_MODES);
	im_k       = dmatrix(Nfreq, MAX_MODES);
	mode_S     = dmatrix(Nfreq, MAX_MODES);
	mode_R     = dmatrix(Nfreq, MAX_MODES);
	read_dispersion_file();

	f_step = f_vec[ 2 ] - f_vec[ 1 ];

  // if the option --nfft was not passed to the main program then the 
  // NFFT default was 0. Otherwise it has the requested positive value. 
  // Here we make sure that whatever the current value of NFFT is 
  // it's not less than 4*n_freqs
  if (NFFT < (4*Nfreq)) {
    NFFT = 4 * Nfreq;
    printf("Minimum NFFT is set to NFFT = %d\n", NFFT);
  }
  
  // if you're given a source file, make sure you can fit the whole thing in NFFT points
  if (source_type.compare("waveform") == 0) {
    int filelines = NCPA::count_rows_arbcol( source_file.c_str() );
    while (filelines > NFFT) {
      NFFT *= 2;
      std::cout << "File " << source_file << " has " << filelines << " lines, increasing NFFT to "
                << NFFT << std::endl;
    }
  }
  
}


NCPA::ModalBroadbandPropagator::~ModalBroadbandPropagator() {
	delete [] mode_count;
	free_dmatrix(re_k,   Nfreq, MAX_MODES);
	free_dmatrix(im_k,   Nfreq, MAX_MODES);								
	free_dmatrix(mode_S, Nfreq, MAX_MODES);
	free_dmatrix(mode_R, Nfreq, MAX_MODES);
}


int NCPA::ModalBroadbandPropagator::read_dispersion_file() {

  // this function reads the file created with write_dispersion()
  // reads freq #modes rho(@z=src), rho(@z=rcv), k_horiz, modes(@z=src), modes(@z=receiverheight)
  int m,n;
  double x1,x2,y1,y2;
  FILE *fp = fopen(dispersion_input_file.c_str(),"r");
  if(fp==NULL){
      std::ostringstream es;
      es << "file << " << dispersion_input_file.c_str() << " could not be opened.\n";
      throw std::invalid_argument(es.str());
  }
  //read_header(f); // use it if file has header
  printf("--> Opening dispersion file %s ...if the program chokes here the file may be corrupt.\n",\
                             dispersion_input_file.c_str());
  // read line by line;
  f_vec[0] = 0; // the DC frequency
  mode_count[0] = 0;
  n=1; // start from 1 because f_vec[0] = 0
  //n = 0;
  while(fscanf(fp,"%lf",&f_vec[n])==1){
      fscanf(fp,"%d",&mode_count[n]);
      fscanf(fp,"%lf",&rho_zsrc);
      fscanf(fp,"%lf",&rho_zrcv);
      for(m=0;m<mode_count[n];m++){
          if(fscanf(fp,"%lf %lf %lf %lf",&x1,&x2,&y1,&y2)!=4){
              fprintf(stderr,"%s: read_dispersion_data format error\n",dispersion_input_file.c_str());
              exit(1);
          }
          re_k[n][m]  =x1;
          im_k[n][m]  =x2;
          mode_S[n][m]=y1;
          mode_R[n][m]=y2;
      }
      ++n;
  }
  fclose(fp);
  printf("--> Found modes at %d positive frequencies in file %s\n",n, dispersion_input_file.c_str());

  return 0; 
}


int NCPA::ModalBroadbandPropagator::calculate_waveform() {
	int i,n;
  double rr, tskip, fmx, t0;	
  std::complex<double> cup, *dft_vec, *pulse_vec, *arg_vec;
  //complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  

  // if zero attenuation requested set the imaginary part of the wavenumber = 0
  if (zero_attn_flag) {
    for (n=0; n<Nfreq; n++) {
      for(int m=0;m<mode_count[n];m++){
          im_k[n][m]  = 0.0;
      }
    }
  }

  dft_vec   = new std::complex<double> [Nfreq];
  transfer_function = new std::complex<double> [Nfreq];
  pulse_vec = new std::complex<double> [NFFT];
  arg_vec   = new std::complex<double> [NFFT]; 

  fmx = ((double)NFFT)*f_step; // max frequency
								      
  get_source_spectrum( dft_vec, pulse_vec, arg_vec );
	
	// if (Nr == 1) { // propagation to one receiver at distance RR from source
	//     std::cout << "--> Doing pulse propagation source-to-receiver at one range: " 
	//          << r_vec[0]/1000.0 << " km" << std::endl;
	         
	//     t0 = r_vec[0] /max_cel;
 //      std::memset( transfer_function,  0, Nfreq * sizeof( std::complex< double > ) );
 //      compute_modal_sum( r_vec[ 0 ] );
	//     //
	//     // fft propagation: note that here 'pulse_vec' is overwritten with the propagated pulse
	//     //							
 //      fft_pulse_prop( t0, r_vec[0], dft_vec, pulse_vec);						      
      
 //      // save propagated pulse to file
	//   // DV 20150514 - factor of two necessary because we only used 
	//   // the positive freq. spectrum in the fft
 //      // DV 20170810 - parameter factor to make it easy to agree with other codes 
 //      // (e.g. Roger Waxler's modal code)
 //      double factor = 2.0;													
 //      f = fopen(waveform_out_file.c_str(),"w");
 //      for(i=0; i<NFFT; i++) {
 //          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx+t0, factor*real(pulse_vec[i]));
 //      }
 //      fclose(f);

 //      std::cout << "--> Propagated pulse saved in file: '" << waveform_out_file << "'" << std::endl
 //           << "' with format: | Time (s) | Re(pulse) |" << std::endl;  

	// }
	// else { // propagation to several receivers
		// @todo combine these, only difference is the use of multiple ranges?

	  std::cout << "--> Propagating pulse from source-to-receivers on grid ..." << std::endl;					      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(waveform_out_file.c_str(),"w");
      tskip = 0.0;
      //double DR = r_vec[1] - r_vec[0];
      //for(n=0; n<=(int)(std::floor((r_vec[Nr-1]-r_vec[0])/DR)); n++) {
      for (n=0; n < Nr; n++) {
          //rr= R_start + DR*n;
      	  rr = r_vec[ n ];
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

          //std::memset( transfer_function,  0, Nfreq * sizeof( std::complex< double > ) );
          std::fill( transfer_function, transfer_function + Nfreq, std::complex<double>{} );
          compute_modal_sum( rr );

          // the call with NFFT as argument			            
          fft_pulse_prop( t0, rr, dft_vec, pulse_vec );

          // DV 20170810 - parameter 'factor' to make it easy to agree with other codes 
          // (e.g. Roger Waxler's modal code)
          double factor = 2.0;								              
          for(i=0;i<NFFT;i++){
              fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx+t0, factor*real(pulse_vec[i]));
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", NFFT, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", waveform_out_file.c_str());
      printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
  // }

  delete [] pulse_vec;
  delete [] arg_vec;
  delete [] dft_vec;
  delete [] transfer_function;

  return 0;
}


void NCPA::ModalBroadbandPropagator::compute_modal_sum( double range ) {
  int i, smooth_space, j;
  std::complex<double> I (0.0, 1.0);
  std::complex< double > k_H, cup, pot, mode_prod;
  double sqrt_rho_ratio = sqrt(rho_zrcv/rho_zsrc);
  std::complex<double> expov8pir = I*exp(-I*Pi*0.25)/sqrt(8.0*Pi*range);

  for(i=0;i<Nfreq;i++){
      cup=0.0;
      for(j=0;j<mode_count[i];j++){
          k_H = re_k[i][j]+I*im_k[i][j];
          pot = exp(I*range*k_H);
          mode_prod = mode_S[i][j]*mode_R[i][j];
          pot = mode_prod*pot;
          pot = pot/sqrt(k_H);
          cup = cup+pot; // modal sum: sum( exp(ikr)/sqrt(k)*Vr*Vs )
      }
      // Note on mode scaling in this code vs. the modes in Computational Oc. Acoust. 1994: 
      // V_book =  sqrt(rho)*V_in_this_code; 
      // Thus the formula for the Fourier component of pressure using the modes in this code
      // is given in DV Modess notes eq. 25 pg. 3 and contains the factor sqrt_rho_ratio
      // as opposed to the factor 1/rho(z_s) in the book eq. 5.14 pg 274.
      transfer_function[i] = expov8pir*cup*sqrt_rho_ratio; // note sqrt_rho_ratio
      //arg_vec[i0+i]=sqrt(I/(8.0*Pi*range))*cup;

  }

  smooth_space=(int)floor(0.1*Nfreq); // smoothly zero out on right; as in RW (July 2012)

  for(i=Nfreq-smooth_space;i<Nfreq;i++){
      transfer_function[i]=transfer_function[i]*half_hann(Nfreq-smooth_space,Nfreq-1,i); // changed df to 1 as df doesn't make sense here
  }
}

