#include "EPadeBroadbandPropagator.h"


namespace NCPA {

	EPadeBroadbandPropagator::EPadeBroadbandPropagator( ParameterSet *param ) {
		// Initialize the frequency parameters
		f_min = param->getFloat( "f_min" );
		f_max = param->getFloat( "f_max" );
		f_step = param->getFloat( "f_step" );
        Nfreq = (int)( std::floor( ( f_max - f_min ) / f_step ) ) + 2;
		NFFT = 4 * Nfreq; 
		// f_vec  = NCPA::zeros<double>( Nfreq );
		f_vec.resize( Nfreq );
        for ( int fi = 1; fi < Nfreq; fi++ ) {
			f_vec[ fi ] = f_min + ( (double)fi - 1) * f_step;
        }
		f_center = f_vec[ Nfreq - 1] / 5; // @todo: Make this an optional parameter 
		// transfer_function = NCPA::zeros<std::complex<double>>( Nfreq );
		transfer_function.resize( Nfreq );
		
		// Other parameters
		source_type = param->getString( "source" );
		std::cout << "Source type: " << source_type << std::endl; 
		receiver_range_km = param->getFloat( "receiver_range_km" );
		rr = receiver_range_km * 1000.0; // Convert to meters
		max_cel = param->getFloat( "max_celerity" );

		t0 = rr / max_cel;

		EPadeSolver *solver = new EPadeSolver( param );
		solver->solve( transfer_function );
		std::cout << "Transfer function calculated." << std::endl;

		// scale the transfer function by 1 / 4*1000*PI to agree with ModBB
		const double _factor = 4000.0 * PI;
		for (int i = 1; i < Nfreq; i++) {
			transfer_function[i] = transfer_function[i] / _factor;
		}

		delete solver;
	}

	EPadeBroadbandPropagator::~EPadeBroadbandPropagator() {
		// @todo: Double check if all allocated memory is freed
	}

	int EPadeBroadbandPropagator::calculate_waveform() {
		// Apply the half-hann window to the transfer function as done in ModBB
		int smooth_space=(int)floor(0.1*Nfreq); // smoothly zero out on right; as in RW (July 2012)


		for(int i=Nfreq-smooth_space;i<Nfreq;i++){
			transfer_function[i]=transfer_function[i]*half_hann(Nfreq-smooth_space,Nfreq-1,i); // changed df to 1 as df doesn't make sense here
		}

		// Create the arrays needed for pulse propagation 
  		// std::complex<double> *dft_vec, *pulse_vec, *arg_vec;
		// dft_vec   = new std::complex<double> [Nfreq];
		// pulse_vec = new std::complex<double> [NFFT];
		// arg_vec   = new std::complex<double> [NFFT];
		std::vector<std::complex<double>> dft_vec( Nfreq );
		std::vector<std::complex<double>> pulse_vec( NFFT );
		std::vector<std::complex<double>> arg_vec( NFFT );


		double fmx = ((double)NFFT)*f_step; // max frequency
		
		std::cout << "Starting get_source_spectrum()" << std::endl;
  		get_source_spectrum( dft_vec, pulse_vec, arg_vec );

		std::cout << "Starting fft_pulse_prop()" << std::endl;
		fft_pulse_prop( t0, rr, dft_vec, pulse_vec );

		 // @todo: Remove this after testing 
		// output the transfer function to file 
		std::cout << "--> Saving transfer function to file: 'transf.pe'" << std::endl;
		FILE *f_transf = fopen("transf.pe","w");
		for (int i = 0; i < Nfreq; i ++) {
			fprintf(f_transf, "%5.6f %15.6e %15.6e\n", f_vec[i], real(transfer_function[i]), imag(transfer_function[i]));
		}
		fclose(f_transf);


		std::cout << "--> Saving propagated pulse to file: 'waveform.pe'" << std::endl;
		FILE *f_pulse = fopen("waveform.pe","w");
		double factor = 2.0; // Same factor as in ModBB, but we also implicitly multiply with I to agree with ModBB convention
		// Note implicit multiplication with I to agree with ModBB convention
		for(int i=0; i<NFFT; i++){
			fprintf(f_pulse,"%10.3f %12.6f %15.6e %15.6e\n", rr/1000.0, 1.0*i/fmx+t0, -factor*imag(pulse_vec[i]), factor*real(pulse_vec[i]));
		}
		fclose(f_pulse);



		// delete[] arg_vec; 
		// delete[] pulse_vec;
		// delete[] dft_vec; 

		return 0;
	}

}