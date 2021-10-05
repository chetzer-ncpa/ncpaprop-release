#ifndef NCPAPROP_BROADBANDPROPAGATOR_H_INCLUDED
#define NCPAPROP_BROADBANDPROPAGATOR_H_INCLUDED

#include "parameterset.h"
#include <complex>

namespace NCPA {

	class BroadbandPropagator {

	public:
		// BroadbandPropagator( ParameterSet *param );
		virtual ~BroadbandPropagator();

		virtual int calculate_waveform() = 0;

		

	protected:
		bool zero_attn_flag, single_receiver;
		int    Nfreq, Nr, NFFT, src_flg; //, *mode_count;
		double f_center, f_step, max_cel; //, R_start, R_end, DR;
		// double rho_zsrc, rho_zrcv;
		double *f_vec, *r_vec;
		// double **re_k, **im_k, **mode_S, **mode_R;
		std::complex< double > *transfer_function;
		std::string waveform_out_file, source_type, source_file;
		// std::string dispersion_input_file;

		// int read_dispersion_file();
		// void compute_modal_sum( double range );
		int get_source_spectrum( std::complex<double> *dft_vec, std::complex<double> *pulse_vec, 
			std::complex<double> *arg_vec );
		std::complex<double> pulse_spec_fit(double scale, double x);
		void model_pulse_fft(double power, double scale, double dt, std::complex<double> *dft_vec);	
		/* Power law pulse model with exponential decay. */
		double model_pulse_shape(double power,double scale,double x);
		int load_source_spectrum(double *freqv, std::complex<double> *dft_vec);
		int load_source_pulse_td(std::vector<double> &t, std::vector<double> &tdp );
		void fft_pulse_prop( double t0, double range, std::complex<double> *dft_vec, 
			std::complex<double> *pulse_vec );
		double half_hann(int begin,int end,int i);

	};
}









#endif
