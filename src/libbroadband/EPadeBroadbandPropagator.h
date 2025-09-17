#ifndef NCPAPROP_EPADEBROADBANDPROPAGATOR_H_INCLUDED
#define NCPAPROP_EPADEBROADBANDPROPAGATOR_H_INCLUDED

#include "parameterset.h"
#include "BroadbandPropagator.h"
#include <complex>

#include "EPadeSolver.h"

namespace NCPA {

	class EPadeBroadbandPropagator : public BroadbandPropagator {

	public:
		EPadeBroadbandPropagator( ParameterSet *param );
		~EPadeBroadbandPropagator();

		// int calculate_waveform( std::complex<double> *transf );
		int calculate_waveform();

	protected:
		// Arrays
		// double *f_vec = nullptr; 
		// std::complex<double> *transf = nullptr;

		// Doubles 
		double f_min, f_max; //, f_step, f_center; are defined in parent
		double receiver_range_km, rr; // rr is the range in meters
		double max_cel; // @todo: Include this as an optional parameter (is set to 340 atm)
		double t0;

		// Ints 
		// int Nfreq; // Defined in parent
		// int NFFT; // Number of FFT points

		// Strings 
		// std::string source_type;


	};

}

#endif