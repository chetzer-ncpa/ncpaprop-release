#ifndef NCPAPROP_MODALBROADBANDPROPAGATOR_H_INCLUDED
#define NCPAPROP_MODALBROADBANDPROPAGATOR_H_INCLUDED

#include "parameterset.h"
#include "BroadbandPropagator.h"
#include <complex>

namespace NCPA {

	class ModalBroadbandPropagator : public BroadbandPropagator {

	public:
		ModalBroadbandPropagator( ParameterSet *param );
		~ModalBroadbandPropagator();

		int calculate_waveform();

	protected:
		int *mode_count;
		double rho_zsrc, rho_zrcv;
		double **re_k, **im_k, **mode_S, **mode_R;
		// std::complex< double > *modal_sum;
		std::string dispersion_input_file;

		int read_dispersion_file();
		void compute_modal_sum( double range );

	};

}



#endif
